using MeasurementBrowser

# Get directory from command line or use default
if length(ARGS) > 0
    measurement_dir = ARGS[end]
else
    measurement_dir = "/home/dgfl/work/Borg/RuO2 processing/v2 (liftoff)/099_MeasData/"
end

if !isdir(measurement_dir)
    println("Error: Directory '$measurement_dir' does not exist!")
    return
end

app = start_browser(measurement_dir)
