using DataBrowserAnnotations
import CImGui as ig

"""Show the tag-loading failure without blocking the rest of the browser."""
function _render_tag_state_error!(state::BrowserState)::Nothing
    message = state.tag_state_error
    isempty(message) && return nothing
    ig.TextColored((1.0, 0.5, 0.5, 1.0), "tags error")
    if ig.BeginItemTooltip()
        ig.PushTextWrapPos(ig.GetFontSize() * 35.0)
        ig.TextUnformatted(message)
        ig.PopTextWrapPos()
        ig.EndTooltip()
    end
end

"""Apply the dominant tag color and report whether the caller must pop it."""
function _push_tag_text_style!(
    tag_state::DataBrowserAnnotations.Tags.TagState,
    key::AbstractString,
    ancestor_keys::Vector{String},
)::Bool
    effective_tags = DataBrowserAnnotations.Tags.effective(tag_state, key, ancestor_keys)
    color = DataBrowserAnnotations.Tags.dominant_color(tag_state, effective_tags)
    color === nothing && return false
    red = Float32(color[1]) / 255f0
    green = Float32(color[2]) / 255f0
    blue = Float32(color[3]) / 255f0
    ig.PushStyleColor(ig.ImGuiCol_Text, (red, green, blue, 1.0f0))
    return true
end
