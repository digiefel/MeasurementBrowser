using DataBrowser
using DelimitedFiles: readdlm
import DataBrowser:
    collection,
    entries,
    fingerprint,
    id,
    item_data,
    label,
    metadata,
    process,
    project_name,
    read,
    reconstruct,
    source_id,
    source_item_path,
    source_items,
    source_label

length(ARGS) == 1 || error("Usage: julia --project project.jl DATA_DIRECTORY")

struct MicrographDirectory <: AbstractDataSource
    root::String
end

struct MicrographFile <: AbstractDataSourceItem
    path::String
    modified::Float64
end

source_id(source::MicrographDirectory)::String = abspath(source.root)
source_label(source::MicrographDirectory)::String = basename(abspath(source.root))

# The engine scans with `cancel_token`, `on_progress`, and `on_item` keywords so a slow source can
# stream and be interrupted. A source that discovers everything at once absorbs and ignores them:
# the engine queues the returned batch when `on_item` was never called.
function source_items(source::MicrographDirectory; kwargs...)::Vector{MicrographFile}
    paths = sort!(filter(
        path -> endswith(lowercase(path), ".txt"),
        readdir(source.root; join=true),
    ))
    return [MicrographFile(path, stat(path).mtime) for path in paths]
end

id(file::MicrographFile)::String = abspath(file.path)
label(file::MicrographFile)::String = basename(file.path)
source_item_path(file::MicrographFile)::String = file.path
fingerprint(file::MicrographFile)::Float64 = file.modified

# A typed project supplies its own collection levels too. `collection` returns `AbstractCollection`
# values, not names: the level is a value the project can dispatch on later, and the package derives
# its identity and label from the contract below.
struct MicrographSet <: AbstractCollection
    name::String
end

id(set::MicrographSet)::String = set.name
label(set::MicrographSet)::String = set.name
# Required for collections: the index stores a level's kind, identity, and metadata, never the
# value, and rebuilds it from those. Items get a slow fallback here; collections do not.
reconstruct(::Type{MicrographSet}, identity::AbstractString, ::Dict)::MicrographSet =
    MicrographSet(String(identity))

struct Micrograph <: AbstractDataItem
    name::String
    pixels::Matrix{Float32}
    exposure_ms::Float64
end

label(image::Micrograph)::String = image.name
collection(::Micrograph)::Vector{AbstractCollection} =
    AbstractCollection[MicrographSet("Micrographs")]
metadata(image::Micrograph)::Dict{Symbol,Any} = Dict{Symbol,Any}(
    :exposure_ms => image.exposure_ms,
    :height_px => size(image.pixels, 1),
    :width_px => size(image.pixels, 2),
)
item_data(image::Micrograph)::Matrix{Float32} = image.pixels

function process(image::Micrograph)::Micrograph
    low, high = extrema(image.pixels)
    scale = high == low ? one(Float32) : high - low
    normalized = (image.pixels .- low) ./ scale
    return Micrograph(image.name, normalized, image.exposure_ms)
end

# `read` is the only stage that touches the source; `entries` is a pure function of its result.
read(::MicrographDirectory, file::MicrographFile)::Matrix{Float32} =
    Float32.(readdlm(file.path, ','))

entries(file::MicrographFile, pixels::Matrix{Float32})::Vector{Micrograph} =
    [Micrograph(splitext(basename(file.path))[1], pixels, 10.0)]

struct MicrographProject <: DataBrowser.AbstractProject end

project_name(::MicrographProject) = "Typed micrographs"

project = MicrographProject()
source = MicrographDirectory(only(ARGS))
workspace = open_workspace(project, source)
open_browser(workspace)
