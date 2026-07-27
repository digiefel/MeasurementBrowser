"""
Supertype of every DataBrowser project.

A project answers the declarations in this file and implements the pipeline stages in
`stage_contract.jl`. It carries no package-owned state: caches, jobs, indexes, and browser state
belong to the workspace. The registration dialect's recipe-holding project is one implementation
of this contract, not the contract itself.
"""
abstract type AbstractProject end

"""Return the stable name used to identify a project."""
function project_name end

"""Return a short human-readable description of a project."""
function project_description end
