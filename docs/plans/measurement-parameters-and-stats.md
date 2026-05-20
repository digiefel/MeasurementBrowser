# Measurement Metadata Contract

## Goal

Make measurement metadata simple enough that users can trust it and developers can extend it without
guessing.

Every value shown by the browser must answer three questions clearly:

```text
What object does this describe?
Where did it come from?
Was it parsed directly or computed later?
```

This document defines the user-facing contract. It does not prescribe the internal implementation.
Here, "user" includes both people using the browser and people adding a new project or measurement
procedure.

## Buckets

Measurement metadata has three public buckets.

```text
Device parameters
  Facts about the physical device.
  Example: area, film thickness, coordinates.

Measurement procedure parameters
  Facts about how this measurement entry was configured or identified.
  Parsed cheaply from the filename, CSV header, or the local file structure.

Analyzed stats
  Facts computed from loaded data, cached analysis, or measurement history.
```

The same concept must not appear in multiple buckets with competing meanings.

## User Contract

- Device parameters describe the device, not the measurement.
- Measurement procedure parameters describe the measurement entry before analysis.
- Analyzed stats describe what the data or measurement history shows after analysis.
- Project/procedure authors must have one obvious place to declare the fields their procedure
  exposes.
- Project/procedure authors must be able to tell whether a field should be parsed as a procedure
  parameter or computed as an analyzed stat without reading UI, cache, or plotting code.
- Missing numeric values are explicit and consistent, not silently absent because one file type took
  a different code path.
- Adding an ordinary field should require one obvious declaration and one obvious test.
- UI and cache should preserve the bucket a value belongs to.
- Tests are written before implementation and define the intended shape.

## First Contract: PUND Family

PUND, wakeup, and fatigue are the first family to migrate because they currently mix filenames,
headers, selected entries, and computed values.

For this family, the browser must expose:

```text
Analyzed stats
  pulse history and actual values computed from the selected waveform
```

The first required history stats for every PUND-family logical measurement are:

```text
wakeup_count
wakeup_f
wakeup_V
fatigue_count
fatigue_f
fatigue_V
```

The count fields describe pulses already applied to the device before the readout waveform,
including pulses sent earlier in the same file. Count fields default to `0` when that kind of pulse
history is absent. Frequency and voltage fields default to `NaN` when they are unavailable.

Wakeup and fatigue files have different structures and must have separate extractors. They share
the public field names, not a generic parsing path. Wakeup counts accumulate per device and wakeup
condition across files. Fatigue files can contain a readout at fatigue count `0`; that is the first
fatigue logical measurement, not "no fatigue file". Fatigue measurements inherit prior wakeup
history for the same device, and later standalone PUND measurements inherit the latest fatigue
history for that device.

The first required analyzed voltage stats for every PUND-family logical measurement are:

```text
V_base
V_min
V_max
V_amp
```

These are stats because they describe the loaded waveform, not merely what the measurement was
requested to do.

The existing fixture tests are the contract for the first migration.

## Non-Goals

- Do not keep compatibility aliases for known-wrong PUND metadata names.
- Do not move device metadata into measurement parsing.
- Do not treat waveform-derived values as parsed measurement parameters.
- Do not make UI text assertions the source of truth for metadata correctness.

## Rollout

1. Add contract tests for the target metadata shape using real fixture files. Done.
2. Add the metadata storage fields to the data model. Done.
3. Simplify the scan path so project parsers create `MeasurementInfo` directly. Done.
4. Keep one full-list post-processing step for history-derived metadata. Done.
5. Move waveform/history-derived values into analyzed stats. Done for RuO2 PUND.
6. Migrate PUND/wakeup/fatigue first and delete old mixed keys. Done in source code.
7. Use the same pattern for other procedure families later.

## Current State

Project parsers now return `MeasurementInfo` records directly. `MeasurementItem` was removed, and
measurement identity is exposed as `unique_id`.

RuO2 PUND-family logical measurements use one full-list post-processing step for wakeup/fatigue
history and waveform stats after all measurements for a device are known.

RuO2 PUND-family history now lives in `MeasurementInfo.stats`: `wakeup_count`, `wakeup_f`,
`wakeup_V`, `fatigue_count`, `fatigue_f`, and `fatigue_V`.

RuO2 PUND-family waveform stats also live in `MeasurementInfo.stats`: `V_base`, `V_min`, `V_max`,
`V_amp`, `frequency_kHz`, and `Pr_max_uCcm2`. Voltage stats are rounded to one decimal place.

PUND-adjacent loading paths use local selection parameters only to find the selected waveform. They
do not reconstruct public wakeup/fatigue history in UI or plotting code.

PUND-adjacent loading paths use the source file kind to decide how to load the waveform. They do
not use `fatigue_count > 0` as a proxy for fatigue files, because count `0` is valid inside a
fatigue file and standalone PUND can still have positive fatigue history.

The cache metadata field list now uses the new PUND fields. Caches written with old PUND field names
are treated as out of date so the normal rebuild path can replace them.

The info panel now renders `MeasurementInfo.parameters` and `MeasurementInfo.stats` directly instead
of hardcoding PUND-specific stat computation in the UI.

RuO2 PUND stats now have one file/dataframe entry point. It computes `V_base`, `V_min`, `V_max`,
`V_amp`, `frequency_kHz`, and `Pr_max_uCcm2` when the required device area is available.

## Open Question

How much of this should align now with the future project/procedure API rewrite?
