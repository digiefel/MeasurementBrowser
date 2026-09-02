# Bench-only source: one real file presented as many items (distinct ids, shared path).

using DataBrowserAPI: AbstractDataSource, AbstractDataSourceItem
import DataBrowserAPI:
    fingerprint, id, label, metadata, source_id, source_item_path, source_items, source_label
using DataBrowserSources: FileFingerprint, file_fingerprint

struct BenchFile <: AbstractDataSourceItem
    filepath::String
    filename::String
    relative_path::String
    fingerprint::FileFingerprint
end

id(file::BenchFile) = file.relative_path
label(file::BenchFile) = file.relative_path
fingerprint(file::BenchFile) = file.fingerprint
source_item_path(file::BenchFile) = file.filepath
metadata(file::BenchFile) = Dict{Symbol,Any}(:filename => file.filename)

struct BenchSource <: AbstractDataSource
    identity::String
    files::Vector{BenchFile}
end

source_id(source::BenchSource) = source.identity
source_label(source::BenchSource) = source.identity
source_items(source::BenchSource; kwargs...) = source.files
Base.copy(source::BenchSource) = BenchSource(source.identity, copy(source.files))

"""`n` items that all read `template`. `relpath(i)` is the fake id."""
function alias_file(template::AbstractString, n::Integer, relpath)::Vector{BenchFile}
    path = abspath(template)
    isfile(path) || error("Missing template $path")
    token = file_fingerprint(path)
    return BenchFile[let rel = String(relpath(i))
        BenchFile(path, basename(rel), rel, token)
    end for i in 1:n]
end
