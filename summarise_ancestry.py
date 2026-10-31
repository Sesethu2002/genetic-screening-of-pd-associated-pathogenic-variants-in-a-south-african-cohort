import csv
import statistics

components = ["NAMA", "EUR", "AFR", "EAS", "MAL"]
data = []

with open("cohort_ancestry.tsv", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if len(row) != 7:
            print("Skipping unexpected row:", row)
            continue

        sample_id, diagnosis = row[0], row[1]

        try:
            values = [float(x) * 100 for x in row[2:7]]
        except ValueError:
            print("Skipping non-numeric row:", row)
            continue

        data.append((diagnosis, values))

groups = {
    "Overall": data,
    "PD": [r for r in data if r[0] == "PD"],
    "Control": [r for r in data if r[0] == "Control"]
}

output = "global_ancestry_summary.tsv"

with open(output, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "Group", "Ancestry", "N", "Mean_percent",
        "Median_percent", "SD_percent",
        "Minimum_percent", "Maximum_percent"
    ])

    for group_name, records in groups.items():
        for i, component in enumerate(components):
            vals = [r[1][i] for r in records]
            writer.writerow([
                group_name,
                component,
                len(vals),
                round(statistics.mean(vals), 4),
                round(statistics.median(vals), 4),
                round(statistics.stdev(vals), 4) if len(vals) > 1 else 0,
                round(min(vals), 4),
                round(max(vals), 4)
            ])

print("Summary saved to", output)
print("Individuals included:")
for name, records in groups.items():
    print(name + ":", len(records))
