"""Optional hook for cooperative scan cancellation (wired by the engine)."""
const SCAN_CANCEL_CHECK = Ref{Function}(() -> nothing)

"""Install the engine's cancellation probe for directory scans."""
function set_scan_cancel_check!(check::Function)::Nothing
    SCAN_CANCEL_CHECK[] = check
    return nothing
end

function _check_scan_cancel()::Nothing
    SCAN_CANCEL_CHECK[]()
    return nothing
end
