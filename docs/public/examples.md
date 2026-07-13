# Examples

The repository's `examples/` directory contains complete projects arranged by the amount of API
they introduce.

1. **Simple CSV** — one file becomes one item; only `read` is needed.
2. **Several items from one source** — `entries` separates loaded data and supplies item metadata.
3. **Typed data** — a concrete `AbstractDataItem` adds domain behavior through multiple dispatch.
4. **Domain project** — several registrations, collections, processing, and analysis work together.

Read the examples in order when learning the API. Each project introduces one additional reason to
write project code; none begins with the full interface.

The examples are also useful as templates. Copy the smallest one matching the structure of your
data, then replace its source interpretation with your own.
