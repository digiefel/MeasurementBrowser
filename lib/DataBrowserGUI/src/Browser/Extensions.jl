"""Dispatch-based GUI extension surface for optional browser panels and menus."""
abstract type GuiExtension end

extension_id(ext::GuiExtension)::String = String(nameof(typeof(ext)))
init!(ext::GuiExtension, state) = nothing
menu!(ext::GuiExtension, state) = nothing
draw!(ext::GuiExtension, state) = nothing
reset!(ext::GuiExtension, state) = nothing
shutdown!(ext::GuiExtension, state) = nothing
is_ready(ext::GuiExtension, state)::Bool = true
save_view(ext::GuiExtension, state)::Dict{String,Any} = Dict{String,Any}()
load_view!(ext::GuiExtension, state, view::Dict{String,Any}) = nothing
open_detached_plot!(ext::GuiExtension, state; kwargs...)::Bool = false

const _GUI_EXTENSION_TYPES = Type{<:GuiExtension}[]

"""Register an extension type; load order follows registration order."""
function register_gui_extension!(::Type{T}) where {T<:GuiExtension}
    filter!(entry -> entry != T, _GUI_EXTENSION_TYPES)
    push!(_GUI_EXTENSION_TYPES, T)
    return nothing
end

"""Instantiate one extension object per registered type."""
function _instantiate_extensions()::Vector{GuiExtension}
    return [T() for T in _GUI_EXTENSION_TYPES]
end
