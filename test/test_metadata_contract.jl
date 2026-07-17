const METADATA_CONTRACT_API = DataBrowserAPI
struct ContractSourceItem <: METADATA_CONTRACT_API.AbstractDataSourceItem
    key::String
end

METADATA_CONTRACT_API.id(item::ContractSourceItem) = item.key
METADATA_CONTRACT_API.label(item::ContractSourceItem) = item.key

struct MetadataContractDataItem <: METADATA_CONTRACT_API.AbstractDataItem
    metadata::Dict{Symbol,Any}
end

METADATA_CONTRACT_API.metadata(item::MetadataContractDataItem) = item.metadata

struct MetadataContractSourceItem <: METADATA_CONTRACT_API.AbstractDataSourceItem
    metadata::Dict{Symbol,Any}
end

METADATA_CONTRACT_API.metadata(item::MetadataContractSourceItem) = item.metadata

@testset "source and item metadata contract" begin
    @test METADATA_CONTRACT_API.metadata(ContractSourceItem("source")) == Dict()
    @test METADATA_CONTRACT_API.metadata(
        MetadataContractDataItem(Dict{Symbol,Any}()),
    ) == Dict()

    source_metadata = Dict{Symbol,Any}(:source => "camera")
    item_metadata = Dict{Symbol,Any}(:exposure => 2.5)

    @test METADATA_CONTRACT_API.metadata(MetadataContractSourceItem(source_metadata)) ===
        source_metadata
    @test METADATA_CONTRACT_API.metadata(MetadataContractDataItem(item_metadata)) === item_metadata

    mktempdir() do directory
        path = joinpath(directory, "measurement.csv")
        write(path, "value\n1\n")
        file = DataBrowserSources.index_source_file(path)
        @test METADATA_CONTRACT_API.metadata(file) == Dict{Symbol,Any}(
            :filename => "measurement.csv",
        )
    end
end
