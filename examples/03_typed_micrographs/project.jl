using DataBrowser
using DelimitedFiles: readdlm
import DataBrowser:
    collection,
    data_items,
    fingerprint,
    item_data,
    item_label,
    metadata,
    process,
    source_id,
    source_item_id,
    source_item_label,
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

function source_items(source::MicrographDirectory)::Vector{MicrographFile}
    paths = sort!(filter(
        path -> endswith(lowercase(path), ".txt"),
        readdir(source.root; join=true),
    ))
    return [MicrographFile(path, stat(path).mtime) for path in paths]
end

source_item_id(file::MicrographFile)::String = abspath(file.path)
source_item_label(file::MicrographFile)::String = basename(file.path)
source_item_path(file::MicrographFile)::String = file.path
fingerprint(file::MicrographFile)::Float64 = file.modified

struct Micrograph <: AbstractDataItem
    name::String
    pixels::Matrix{Float32}
    exposure_ms::Float64
end

item_label(image::Micrograph)::String = image.name
collection(::Micrograph)::Vector{String} = ["Micrographs"]
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

function data_items(
    ::Project,
    ::MicrographDirectory,
    file::MicrographFile,
)::Vector{Micrograph}
    pixels = Float32.(readdlm(file.path, ','))
    name = splitext(basename(file.path))[1]
    return [Micrograph(name, pixels, 10.0)]
end

project = define_project("Typed micrographs"; description="Matrix-valued microscopy images")
source = MicrographDirectory(only(ARGS))
workspace = open_workspace(project, source)
open_browser(workspace)
