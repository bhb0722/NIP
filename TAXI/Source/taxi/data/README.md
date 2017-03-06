We downloaded the data in raw form from [Chris Whong's website](chriswhong.com/open-data/foil_nyc_taxi/). A
sample of this raw-form data is shown in raw_taxi_data_sample.csv. (The entirety of the raw data is too large to
upload, so we only include 50 rows of the data set here.) We then aggregated 4 months worth of the raw data into the
file pickups_aggregated_manhattan.csv using the sql scripts ../sql/load_raw_from_csv.sql and
../sql/aggregate_pickups.sql.
