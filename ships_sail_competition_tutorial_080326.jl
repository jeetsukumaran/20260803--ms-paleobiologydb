# Ships that pass in the night: Competition in the fossil record (using Julia)

using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

# Part 1: Get to know the Paleobiology Database (PBDB)!

# Before using the PBDB Julia interface, you will need to check out the actual PBDB website and get a sense of
# how the data are structured. The PBDB is a relational database, which means that data are stored in different tables
# that are linked together by unique identifiers. For example, occurrences are stored in one table,
# while taxonomic information is stored in another table.
# To get a sense of how these tables are linked together, you can use the PBDB website to explore the data and see
# how different tables are related.

# In one of the questions, you are asked to go through your birthplace and look for fossils nearest to your location (the circles)
# Check out both the "General" and "Occurrences" tabs in the dataset, then answer the questions on your activity PDF.

# Using "pbdb_collection", you can pull the exact collection number then call the data to show columns for 
# both location and stratigraphic details:

collection_pk = pbdb_collection("col:178006", show = ["loc", "stratext"], extids = true)

show(describe(collection_pk), allrows = true)

# For the "Occurrences" tab, you can acquire the data using the "pbdb_occurrences" function, as follows:
animal = "Puffinus kanakoffi"

occs = pbdb_occurrences(
    ;
    base_name = animal,
    show      = "full",
    vocab     = "pbdb",
    extids    = true
) # Using Puffinus kanakoffi as an example species, from a location near San Diego
# (11, 134) DataFrame

foreach(println, names(occs))
# There are roughly 134 columns in the occurrences table for P. kanakoffi

describe(occs)
# (134, 7)

show(describe(occs), allrows = true)

combine(groupby(occs, :accepted_rank), nrow)

# Based on our single datum, there are only 11 records within the entire collection, which greatly contrasts the
# collections typically seen in San Diego (our example) with other more diverse collections from around the world.

# Hence, we will call data from a collection near Monterrey, Mexico instead, using data from La Hausteca Canyon

coll_monterrey = pbdb_collection("col:220702", show = "full", extids = true)

show(describe(coll_monterrey), allrows = true)

# Using "pbdb_occurrences", we will call our animal (or phylum) "Mollusca" and our period "Late Aptian"

animal = "Mollusca"
period = "Late Aptian"

occs_mon = pbdb_occurrences(
   base_name = animal,
   interval = period,
   show = "full",
   vocab = "pbdb",
   extids = true
)

foreach(println, names(occs_mon))

describe(occs_mon)

show(describe(occs_mon), allrows = true)

combine(groupby(occs_mon, :accepted_rank), nrow)

clean_taxonomy_flt = row -> row.accepted_rank == "species" || row.accepted_rank == "subspecies"
occs_mon_species = filter(clean_taxonomy_flt, occs_mon) # This is used to filter out the DataFrame to only include occurrences
# for species and sub-species

occs_mon_species[:, r".*_ma.*"] # This is for chronological precision

# Example code to use (You don't have to use all four!)

# println("Direct MA measurements: ", nrow(dropmissing(occs_species, r".*direct_ma_value.*")))
# println("Max MA measurements:    ", nrow(dropmissing(occs_species, r".*max_ma.*")))
# println("Min MA measurements:    ", nrow(dropmissing(occs_species, r".*min_ma.*")))
# println("MA error measurements:  ", nrow(dropmissing(occs_species, r".*ma_error*")))

println("Max MA measurements:    ", nrow(dropmissing(occs_mon_species, r".*max_ma.*")))
println("Min MA measurements:    ", nrow(dropmissing(occs_mon_species, r".*min_ma.*")))

# Only drop records from ":max_ma" and ":min_ma"

occs_ages = dropmissing(occs_species, [:max_ma, :min_ma])

modern_coords = nrow(dropmissing(occs_ages, [:lng, :lat]))
println("Records with modern coordinates:    ", modern_coords)

# These next few pieces of code are mainly used with larger DataFrames for more aggressive filtering (e.g., Carnivora)

# paleo_coords = nrow(dropmissing(occs_ages, [:paleolng, :paleolat]))
# println("Records with paleo-coordinates:     ", paleo_coords)

# complete_spatial = nrow(dropmissing(occs_ages, [:lng, :lat, :paleolng, :paleolat]))
# println("Records with complete spatial data: ", complete_spatial)

# occs_clean = dropmissing(occs_ages, [:lng, :lat, :paleolng, :paleolat])

occs_clean = dropmissing(occs_ages, [:lng, :lat])

species_counts = combine(groupby(occs_clean, :accepted_name), nrow)
genus_counts   = combine(groupby(occs_clean, :genus), nrow)
family_counts  = combine(groupby(occs_clean, :family), nrow)

println("Species:  ", nrow(species_counts))
println("Genera:   ", nrow(genus_counts))
println("Families: ", nrow(family_counts))

# Data filtering should be relaxed for smaller DataFrames (for smaller batches of data)

# Make sure to answer all questions on your activity PDF as you go along the tutorial

# ----------------------------------------------------------------

# Part 2: Exploring different phylogenetic relationships

# You'll need to locate multiple locations, including Poughkeepsie, Saugerties, North of Kingston (Onondaga Formation), and
# Albemarle Sound (Pleistocene)

# When comparing the Kingston fossils to the Albemarle Sound fossils, make sure to write down the names of all phyla/classes present;
# Then, using the PBDB Julia interface, you can query the database for occurrences of these phyla/classes and see
# which phyla/classes have the most species present

# For example, you can use the following code to query the PBDB for occurrences of a general phylum in Onondaga Formation, like Mollusca:

# Using "pbdb_collections" first...

coll_onondaga_1 = pbdb_collection("col:129948", show = "full", extids = true)
show(describe(coll_onondaga_1), allrows = true)

coll_onondaga_2 = pbdb_collection("col:130152", show = "full", extids = true)
show(describe(coll_onondaga_2), allrows = true)

coll_onondaga_3 = pbdb_collection("col:130260", show = "full", extids = true)
show(describe(coll_onondaga_3), allrows = true)

# Then "pbdb_occurrences"

animal = "Mollusca"
period = "Eifelian"

occs = pbdb_occurrences(
   base_name = animal,
   interval = period,
   show = "full",
   vocab = "pbdb",
   extids = true
)

size(occs)
# (1370, 134)

# You can use the "combine" function to see how the records are distributed across different taxonomic ranks:
combine(groupby(occs, :accepted_rank), nrow)
# Row │ accepted_rank  nrow  
#     │ String15       Int64 
#─────┼──────────────────────
#   1 │ genus            385
#   2 │ superfamily       12
#   3 │ subgenus          34
#   4 │ order              7
#   5 │ subclass           9
#   6 │ class             44
#   7 │ species          859
#   8 │ family            16
#   9 │ subspecies         4

clean_taxonomy = row -> row.accepted_rank == "species" || row.accepted_rank == "subspecies"
occs_species = filter(clean_taxonomy, occs)
# (863, 134)

occs_species[:, r".*_ma.*"]
# (863, 14)

println("Max MA measurements: ", nrow(dropmissing(occs_species, r".*max_ma.*")))
println("Min MA measurements: ", nrow(dropmissing(occs_species, r".*min_ma.*")))
# Max MA measurements: 0
# Min MA measurements: 0

# occs_with_ages = dropmissing(occs_species, [:direct_ma_value, :max_ma, :min_ma, :direct_ma_error]) (General code)
occs_with_ages = dropmissing(occs_species, [:max_ma, :min_ma])
# (863, 134)

modern_coords = nrow(dropmissing(occs_with_ages, [:lng, :lat]))
println("Records with modern coordinates:    ", modern_coords)
# Records with modern coordinates:    863

occs_clean = dropmissing(occs_with_ages, [:lng, :lat])
# (863, 134)

# Finally, with the cleaned data, you can use this code to check the number of species, genera, and families present in the dataset:
species_counts = combine(groupby(occs_clean, :accepted_name), nrow)
genus_counts   = combine(groupby(occs_clean, :genus), nrow)
family_counts  = combine(groupby(occs_clean, :family), nrow)

println("Species:  ", nrow(species_counts))
println("Genera:   ", nrow(genus_counts))
println("Families: ", nrow(family_counts))

# Species:  308
# Genera:   162
# Families: 66

# Therefore, there are 308 speices present in the Mollusca dataset set to the Eifelian period, as well as 162 genera and 66 families

# ----------------------------------------------------------------

# Part 3: Creating a diversity curve through time with your phylogenetic relationships

# For the purpose of this tutorial, we will be creating a diversity curve through time for the phylum/class fossil groups
# Brachiopoda, Mollusca, and Bivalvia. We will be using the same data cleaning and processing steps as in others tutorials,
# but we will be creating a diversity loop that takes into account the phylogenetic relationships between these groups, and to
# account for competition between these groups throughout geographic time.

using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

occs_competition_raw_1 = pbdb_occurrences(
    base_name = "Brachiopoda",
    show = "full",
    vocab = "pbdb"
)
# The original DataFrame is (187095 rows × 135 columns), but we will be cleaning our data as such:

occs_competition_1 = dropmissing(occs_competition_raw_1, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

occs_competition_1 = Taxonomy.drop_unqualified_taxa(occs_competition_1, "genus")
# The initial result is (181627 rows x 135 columns)

# Now, to calculate a midpoint age (Ma):

occs_competition_1.midpoint_age = (
    occs_competition_1.max_ma .+ occs_competition_1.min_ma
) ./ 2

interval_competition_1 = combine(
    groupby(occs_competition_1, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)
# 540x3 DataFrame

interval_competition_1.interval_mid_ma = (
    interval_competition_1.interval_max_ma .+ interval_competition_1.interval_min_ma
) ./ 2

interval_competition_1 = sort(interval_competition_1, :interval_mid_ma, rev=true)
# 540x4 DataFrame, including the "interval_mid_ma" column

genera_competition_1 = combine(
    groupby(occs_competition_1, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)
# This DataFrame solely consists of genera belonging to Brachiopoda, and includes the first and last appearance of each genus
# in the fossil record. (3961, 3)

# To create a diversity loop for the diversity curve plot:

diversity = Int[]

for interval_row in eachrow(interval_competition_1)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genera_competition_1)
    )
    push!(diversity, n)
end

interval_competition_1.diversity = diversity

# This generates a 540-element vector via {Int64}; For Mollusca...

using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

occs_competition_raw_2 = pbdb_occurrences(
    base_name = "Mollusca",
    show = "full",
    vocab = "pbdb"
)

# (500659, 135) DataFrame

occs_competition_2 = dropmissing(occs_competition_raw_2, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

occs_competition_2 = Taxonomy.drop_unqualified_taxa(occs_competition_2, "genus")

occs_competition_2.midpoint_age = (
    occs_competition_2.max_ma .+ occs_competition_2.min_ma
) ./ 2

interval_competition_2 = combine(
    groupby(occs_competition_2, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)
# (609, 3) DataFrame

interval_competition_2.interval_mid_ma = (
    interval_competition_2.interval_max_ma .+ interval_competition_2.interval_min_ma
) ./ 2

interval_competition_2 = sort(interval_competition_2, :interval_mid_ma, rev=true)
# (609, 4) DataFrame

genera_competition_2 = combine(
    groupby(occs_competition_2, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)
# (13632, 3) DataFrame

diversity = Int[]

for interval_row in eachrow(interval_competition_2)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genera_competition_2)
    )
    push!(diversity, n)
end

interval_competition_2.diversity = diversity

# For Bivalvia...

using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

occs_competition_raw_3 = pbdb_occurrences(
    base_name = "Bivalvia",
    show = "full",
    vocab = "pbdb"
)

# (210580, 135) DataFrame

occs_competition_3 = dropmissing(occs_competition_raw_3, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

occs_competition_3 = Taxonomy.drop_unqualified_taxa(occs_competition_3, "genus")

occs_competition_3.midpoint_age = (
    occs_competition_3.max_ma .+ occs_competition_3.min_ma
) ./ 2

interval_competition_3 = combine(
    groupby(occs_competition_3, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

# (474, 3) DataFrame

interval_competition_3.interval_mid_ma = (
    interval_competition_3.interval_max_ma .+ interval_competition_3.interval_min_ma
) ./ 2

interval_competition_3 = sort(interval_competition_3, :interval_mid_ma, rev=true)
# (474, 4) DataFrame

genera_competition_3 = combine(
    groupby(occs_competition_3, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)
# (3612, 3) DataFrame

diversity = Int[]

for interval_row in eachrow(interval_competition_3)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genera_competition_3)
    )
    push!(diversity, n)
end

interval_competition_3.diversity = diversity

# We will need to create a separate DataFrame that combines the interval diversity data for all three groups,
# so that we can plot them together on the same graph. We can do this using this code in DataFrames, as follows:

timescale_1 = DataFrame(
    brachiopoda_diversity = interval_competition_1.diversity,
    interval_mid_ma_brach = interval_competition_1.interval_mid_ma
)

timescale_2 = DataFrame(
    mollusca_diversity = interval_competition_2.diversity,
    interval_mid_ma_mollus = interval_competition_2.interval_mid_ma
)

timescale_3 = DataFrame(
    bivalvia_diversity = interval_competition_3.diversity,
    interval_mid_ma_bival = interval_competition_3.interval_mid_ma
)

# Now we can create our diversity curve, as follows:

using Makie, CairoMakie, ColorSchemes

fig = Figure(size = (1750, 550), fontsize = 11.5)
ax1 = Axis(
    fig[1,1],
    xlabel = "Age (Ma)",
    ylabel = "Number of genera",
    title = "Genus diversity of fossils from Brachiopoda, Mollusca, and Bivalvia through time",
    xticks = (interval_competition_1.interval_mid_ma, interval_competition_1.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 4.5,
    xreversed = true,
)
ax2 = Axis(
    fig[1,1],
    xticks = (interval_competition_2.interval_mid_ma, interval_competition_2.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 4.5,
    xreversed = true,
)
ax3 = Axis(
    fig[1,1],
    xticks = (interval_competition_3.interval_mid_ma, interval_competition_3.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 4.5,
    xreversed = true,
)
hidedecorations!(ax2)
hidedecorations!(ax3)

lines!(ax1, interval_competition_1.interval_mid_ma, interval_competition_1.diversity, linewidth = 2, color = :orange)
scatter!(ax1, interval_competition_1.interval_mid_ma, interval_competition_1.diversity, markersize = 8, color = :orange, label = "Brachiopoda")

lines!(ax2, interval_competition_2.interval_mid_ma, interval_competition_2.diversity, linewidth = 2, color = :seagreen3)
scatter!(ax2, interval_competition_2.interval_mid_ma, interval_competition_2.diversity, markersize = 8, color = :seagreen3, label = "Mollusca")

lines!(ax3, interval_competition_3.interval_mid_ma, interval_competition_3.diversity, linewidth = 2, color = :royalblue)
scatter!(ax3, interval_competition_3.interval_mid_ma, interval_competition_3.diversity, markersize = 8, color = :royalblue, label = "Bivalvia")

axislegend(ax1, position = (0.01, 0.925), framevisible = true, labelsize = 11)
axislegend(ax2, position = (0.01, 0.805), framevisible = true, labelsize = 11)
axislegend(ax3, position = (0.01, 0.675), framevisible = true, labelsize = 11)

fig

save("sea_floor_organisms_diversity_curve_080326.png", fig)

# Save your map for reference, and answer the questions listed on your activity PDF