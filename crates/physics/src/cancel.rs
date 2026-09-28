//! Cooperative cancellation of long, non-incremental computations (the factorizations of
//! metal and electrode systems): the game's physics thread runs them inside
//! `cancellable`, and when a newer setup has arrived (another level, an edit) the next
//! `checkpoint` abandons the work by unwinding. A half-built system never reaches a cache:
//! caches insert only complete results, after the computation. Outside `cancellable`
//! (tests, the generator) checkpoints do nothing.

use std::cell::RefCell;

type Check = Box<dyn Fn() -> bool>;

thread_local! {
    static KEEP_GOING: RefCell<Option<Check>> = const { RefCell::new(None) };
}

/// The unwinding payload of a cancelled computation.
struct Cancelled;

/// Abandons the current `cancellable` computation if it is no longer wanted.
pub fn checkpoint() {
    let stop = KEEP_GOING.with(|k| k.borrow().as_ref().is_some_and(|f| !f()));
    if stop {
        // `resume_unwind` does not call the panic hook: no message is printed.
        std::panic::resume_unwind(Box::new(Cancelled));
    }
}

/// Runs `f`, abandoning it (`None`) at a checkpoint where `keep_going` says no. Other
/// panics propagate unchanged.
pub fn cancellable<T>(keep_going: impl Fn() -> bool + 'static, f: impl FnOnce() -> T) -> Option<T> {
    let previous = KEEP_GOING.with(|k| k.replace(Some(Box::new(keep_going))));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f));
    KEEP_GOING.with(|k| *k.borrow_mut() = previous);
    match result {
        Ok(v) => Some(v),
        Err(e) if e.is::<Cancelled>() => None,
        Err(e) => std::panic::resume_unwind(e),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;
    use std::rc::Rc;

    #[test]
    fn stops_at_a_checkpoint_and_passes_results_through() {
        let calls = Rc::new(Cell::new(0));
        let c = calls.clone();
        let r = cancellable(
            move || {
                c.set(c.get() + 1);
                c.get() < 3
            },
            || {
                for _ in 0..10 {
                    checkpoint();
                }
                7
            },
        );
        assert_eq!(r, None);
        assert_eq!(calls.get(), 3);
        assert_eq!(cancellable(|| true, || 7), Some(7));
        // Outside `cancellable` a checkpoint does nothing.
        checkpoint();
    }
}
