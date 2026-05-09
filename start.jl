using MeasurementBrowser

# Get directory from command line or use default
if length(ARGS) > 0
    measurement_dir = ARGS[end]
else
    measurement_dir = "/home/dgfl/work/Borg/RuO2 processing/v2 (liftoff)/099_MeasData/"
end

app = start_browser(measurement_dir)
