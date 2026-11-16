//! Controlled Rayon thread-pool dispatch for host language adapters.

use crate::CoreError;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::sync::{Arc, Mutex};

static CACHED_POOL: Mutex<Option<(usize, Arc<ThreadPool>)>> = Mutex::new(None);

/// Execute `f` using the requested thread budget.
///
/// Setting `nthreads = Some(1)` executes the compiled Rust routine sequentially
/// inside a single-thread pool without touching the global multi-core pool or
/// falling back to interpreted host code.
pub fn run_with_threads<R, F>(nthreads: Option<usize>, f: F) -> Result<R, CoreError>
where
    R: Send,
    F: FnOnce() -> R + Send,
{
    let Some(count) = nthreads else {
        return Ok(f());
    };
    if count == 0 {
        return Err(CoreError("thread count must be at least 1"));
    }
    let pool = {
        let mut guard = CACHED_POOL
            .lock()
            .map_err(|_| CoreError("thread pool mutex was poisoned"))?;
        if let Some((cached_count, ref pool)) = *guard {
            if cached_count == count {
                Arc::clone(pool)
            } else {
                let built = Arc::new(
                    ThreadPoolBuilder::new()
                        .num_threads(count)
                        .build()
                        .map_err(|_| CoreError("failed to initialize Rayon thread pool"))?,
                );
                *guard = Some((count, Arc::clone(&built)));
                built
            }
        } else {
            let built = Arc::new(
                ThreadPoolBuilder::new()
                    .num_threads(count)
                    .build()
                    .map_err(|_| CoreError("failed to initialize Rayon thread pool"))?,
            );
            *guard = Some((count, Arc::clone(&built)));
            built
        }
    };
    Ok(pool.install(f))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_thread_budget_executes_in_rust_with_one_worker() {
        let seen = run_with_threads(Some(1), rayon::current_num_threads).unwrap();
        assert_eq!(seen, 1);
        let seen_two = run_with_threads(Some(2), rayon::current_num_threads).unwrap();
        assert_eq!(seen_two, 2);
    }
}
