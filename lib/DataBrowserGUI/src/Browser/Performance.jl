import ModernGL as gl

"""Read identifying strings from the active OpenGL context."""
function _gl_info()::Dict{Symbol,String}
    try
        return Dict(
            :vendor => unsafe_string(gl.glGetString(gl.GL_VENDOR)),
            :renderer => unsafe_string(gl.glGetString(gl.GL_RENDERER)),
            :version => unsafe_string(gl.glGetString(gl.GL_VERSION)),
            :sl => unsafe_string(gl.glGetString(gl.GL_SHADING_LANGUAGE_VERSION)),
        )
    catch err
        @warn "GL info query failed" error = err
        return Dict{Symbol,String}()
    end
end
