use futures::stream::{AbortHandle, Abortable};
use gloo_worker::HandlerId;
use leptos::prelude::{ArcRwSignal, ArcSignal, Get, Set, StoredValue};
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::{
    cell::RefCell,
    cmp::Ordering,
    future::Future,
    io,
    ops::{Deref, DerefMut},
    sync::LazyLock,
    thread::LocalKey,
    time::Duration,
};
use wasm_bindgen::{prelude::Closure, JsCast};
use web_sys::{self, MediaQueryListEvent};

use utile::{
    drop::ExecuteOnDrop,
    resource::{Compression, RawResource, RawResourceExt, UrlResource},
};

pub static DARK_MODE: LazyLock<ArcSignal<bool>> =
    LazyLock::new(|| media_query_signal("(prefers-color-scheme: dark)"));

pub static PLOTLY_THEME: LazyLock<ArcSignal<&'static plotly::layout::Template>> =
    LazyLock::new(|| {
        ArcSignal::derive(|| {
            if DARK_MODE.get() {
                &*plotly::layout::themes::PLOTLY_DARK
            } else {
                &*plotly::layout::themes::PLOTLY_WHITE
            }
        })
    });

fn media_query_signal(query: &str) -> ArcSignal<bool> {
    let signal = ArcRwSignal::new(false);
    let media_query_handle = on_media_query(query, {
        let signal = signal.clone();
        move |m| signal.set(m)
    });
    StoredValue::new_local(media_query_handle);
    signal.into()
}

fn on_media_query(query: &str, mut f: impl FnMut(bool) + 'static) -> ExecuteOnDrop<impl FnOnce()> {
    let media_query_list = web_sys::window()
        .unwrap()
        .match_media(query)
        .unwrap()
        .unwrap();

    f(media_query_list.matches());

    let f = {
        let media_query_list = media_query_list.clone();
        move |_| f(media_query_list.matches())
    };
    let f: Closure<dyn FnMut(MediaQueryListEvent)> = Closure::wrap(Box::new(f));

    _ = media_query_list.add_event_listener_with_callback("change", f.as_ref().unchecked_ref());

    ExecuteOnDrop::new(move || {
        _ = media_query_list
            .remove_event_listener_with_callback("change", f.as_ref().unchecked_ref());
    })
}

// Horrible hack
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AsJson<T>(pub T);
impl<T> AsJson<T> {
    pub fn inner(self) -> T {
        self.0
    }
}
impl<T> Deref for AsJson<T> {
    type Target = T;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl<T> DerefMut for AsJson<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}
impl<T> Serialize for AsJson<T>
where
    T: Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serde_json::to_string(&self.0)
            .map_err(serde::ser::Error::custom)?
            .serialize(serializer)
    }
}
impl<'de, T> Deserialize<'de> for AsJson<T>
where
    T: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s: &str = serde::Deserialize::deserialize(deserializer)?;
        let v = serde_json::from_str(s).map_err(serde::de::Error::custom)?;
        Ok(Self(v))
    }
}

/// Hack that loads the (compressed!) raw data asyncronously, then deserializes it in-place.
/// Avoids very large allocation of uncompressed data.
// TODO: move resource to uitle.
pub async fn load_large_json<T: DeserializeOwned>(url: &str) -> io::Result<T> {
    struct InMemoryResource {
        data: RefCell<Option<Vec<u8>>>,
    }
    impl RawResource for InMemoryResource {
        const NAMESPACE: &'static str = "memory";
        fn key(&self) -> String {
            unreachable!()
        }
        fn compression(&self) -> Option<Compression> {
            None
        }

        type Reader = std::io::Cursor<Vec<u8>>;
        fn size(&self) -> std::io::Result<u64> {
            unreachable!()
        }
        fn read(&self) -> std::io::Result<Self::Reader> {
            Ok(std::io::Cursor::new(self.data.borrow_mut().take().unwrap()))
        }

        type AsyncReader = std::io::Cursor<Vec<u8>>;
        async fn size_async(&self) -> std::io::Result<u64> {
            unreachable!()
        }
        async fn read_async(&self) -> std::io::Result<Self::AsyncReader> {
            unreachable!()
        }
    }

    let raw = UrlResource::new(url)?.read_vec_async().await?;

    InMemoryResource {
        data: RefCell::new(Some(raw)),
    }
    .decompressed_with(Compression::Brotli)
    .read_json()
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct HandlerIdOrd(HandlerId);
impl std::fmt::Debug for HandlerIdOrd {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <HandlerId as std::fmt::Debug>::fmt(&self.0, f)
    }
}
impl HandlerIdOrd {
    pub fn raw(self) -> usize {
        format!("{self:?}")
            .strip_prefix("HandlerId(")
            .unwrap()
            .strip_suffix(")")
            .unwrap()
            .parse()
            .unwrap()
    }
}
impl PartialOrd for HandlerIdOrd {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for HandlerIdOrd {
    fn cmp(&self, other: &Self) -> Ordering {
        self.raw().cmp(&other.raw())
    }
}

pub async fn wait_for<T>(id: HandlerId, f: impl Fn() -> Option<T>) -> T {
    loop {
        if let Some(t) = f() {
            return t;
        }

        log::info!("[Worker][{id:?}] Data not ready.");
        leptos_ext::util::sleep(Duration::from_millis(200)).await;
    }
}

/// Gives a chance for new requests to replace the current one.
pub async fn yield_now(id: HandlerId, i: Option<usize>) {
    log::info!("[Worker][{id:?}] Yielding: {i:?}");
    utile::time::sleep(Duration::from_millis(0)).await;
}

pub fn spawn_task(
    task: &'static LocalKey<RefCell<Option<(HandlerIdOrd, AbortHandle)>>>,
    id: HandlerId,
    f: impl Future<Output = ()> + 'static,
) {
    let id = HandlerIdOrd(id);
    let (abort_handle, abort_registration) = AbortHandle::new_pair();
    let f = Abortable::new(f, abort_registration);
    let f = async move {
        let _ = f.await;
    };

    task.with(move |task| {
        let task = &mut *task.borrow_mut();
        *task = Some(match task {
            Some((current_id, current_abort_handle)) if *current_id < id => {
                log::info!("[Worker][{id:?}] Aborting {current_id:?}.");
                current_abort_handle.abort();
                wasm_bindgen_futures::spawn_local(f);
                (id, abort_handle)
            }
            Some(_) => {
                log::info!("[Worker][{id:?}] Stale request, dropping.");
                return;
            }
            None => {
                log::info!("[Worker][{id:?}] First request.");
                wasm_bindgen_futures::spawn_local(f);
                (id, abort_handle)
            }
        });
    });
}
