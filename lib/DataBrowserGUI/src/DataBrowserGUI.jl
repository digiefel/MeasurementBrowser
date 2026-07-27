"""CImGui browser shell and panels."""
module DataBrowserGUI

include("Browser.jl")

using .Browser: open_browser, close_browser!, BrowserSession, gui_timings, reset_timings!
export open_browser, close_browser!, BrowserSession, Browser, gui_timings, reset_timings!

end
