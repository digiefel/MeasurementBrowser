# Workspace actions operate on shared data state without depending on the GUI. Session actions
# extend the same generic function, call its workspace method, and synchronize the browser view.
# Keep view updates in this package and workspace behavior in Core. GUI-only actions take a session.
# Both the REPL and GUI can act on the BrowserSession returned by open_browser.

import DataBrowserCore.Workspace: select_items!

"""
    select_items!(session::BrowserSession, items)

Select indexed items by ID, record or data item and reveal their collections in the browser.
Expand their parent collections and scroll to the first selected item. An empty selection clears
the selected items while keeping the current collection view. Tag visibility filters still apply.

This session action calls `select_items!(workspace, items)` and updates the browser view under the
render loop's workspace lock. The workspace action alone selects data without revealing it.
"""
function select_items!(session::BrowserSession, items::AbstractVector)::Nothing
    workspace = session.state.workspace::Workspace.Workspace
    lock(workspace.lifecycle_lock) do
        Workspace.select_items!(workspace, items)
        _reveal_selected_items!(session.state)
    end
    return nothing
end

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
