import csv
import statistics
import math

components = ["NAMA", "EUR", "AFR", "EAS", "MAL"]
data = []

def percentile(values, p):
    values = sorted(values)
    position = (len(values) - 1) * p
    lower = math.floor(position)
    upper = math.ceil(position)

    if lower == upper:
        return values[lower]

    return values[lower] + (values[upper] - values[lower]) * (position - lower)

with open("cohort_ancestry.tsv", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if len(row) != 7:
            continue

        diagnosis = row[1]

        try:
            values = [float(x) * 100 for x in row[2:7]]
        except ValueError:
            continue

        data.append((diagnosis, values))

groups = {
    "Overall": data,
    "PD": [r for r in data if r[0] == "PD"],
    "Control": [r for r in data if r[0] == "Control"]
}

with open("global_ancestry_median_iqr.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "Group", "Ancestry", "N", "Median_percent",
        "Q1_percent", "Q3_percent", "IQR_percent"
    ])

    for group_name, records in groups.items():
        for i, component in enumerate(components):
            vals = [r[1][i] for r in records]

            q1 = percentile(vals, 0.25)
            median = statistics.median(vals)
            q3 = percentile(vals, 0.75)

            writer.writerow([
                group_name,
                component,
                len(vals),
                round(median, 4),
                round(q1, 4),
                round(q3, 4),
                round(q3 - q1, 4)
            ])

print("Summary saved to global_ancestry_median_iqr.tsv")
for name, records in groups.items():
    print(name + ":", len(records))
