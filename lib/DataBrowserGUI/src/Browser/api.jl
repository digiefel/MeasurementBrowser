# Public REPL API for a running browser GUI. Each function takes the `BrowserSession` returned by
# `open_browser`, so the REPL and the GUI act on the same session. (These are the first of these;
# further GUI-facing API belongs here.)
#
# There is one render loop per process, so the always-on timer is a single global (`MAIN_TIMER`);
# the session identifies the GUI and carries the deferred-reset flag.

"""
The always-on main-task timing tree recorded by `@timed` during the render loop. Returns the live
`TimerOutput`, so every TimerOutputs function applies to it (`print_timer`, `flatten`, Tables export,
…). It is lock-free (written only by the render task), so read it while the GUI is idle to avoid
racing a section update.
"""
gui_timings(session::BrowserSession) = MAIN_TIMER

"""
Clear the main-task timer. The render loop performs the reset at the top of its next frame, where no
`@timed` section is open; a no-op if the GUI is not running.
"""
function reset_timings!(session::BrowserSession)::Nothing
    session.state.performance.reset_main_timer = true
    return nothing
end

"""
    wait_browser_ready(session; timeout_s=120) -> session

Wait until the browser has presented a full frame after extension initialization. The preparation
screen does not count. This waits for the GUI, not for all background workspace processing.
Propagate render-task failures; throw if the browser closes or the deadline expires before readiness.
"""
function wait_browser_ready(session::BrowserSession; timeout_s::Real=120)
    deadline = time() + timeout_s
    while !session.state.performance.ready
        if istaskdone(session.task)
            fetch(session.task)
            error("Browser closed before becoming ready")
        end
        time() < deadline || error("Browser did not become ready within $(timeout_s) seconds")
        sleep(0.01)
    end
    return session
end
