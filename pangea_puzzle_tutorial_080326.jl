# The Pangea Puzzle (primary tutorial)

# The main goal of this activity is to track the occurrence of fossil organisms on PaleobiologyDB by creating maps
# and analyzing the maps of each fossil species using Julia coding language

# ------------------------------------------

# Part 1: Construct a map of fossil distributions on the modern continents using Lystrosaurus.

# Set up your environment by opening the following packages
using PaleobiologyDB, PaleobiologyDB.Taxonomy
using DataFrames, Makie, CairoMakie, GeoMakie, ColorSchemes

# Acquire and clean up data for Lystrosaurus fossils

animal = "Lystrosaurus"

lystro_raw = pbdb_occurrences(
    base_name = animal,
    show = ["coords", "classext"],
    vocab = "pbdb",
    extids = true
)

# 124x29 DataFrame; 20 columns and 118 rows omitted

# Print out how many fossil occurrences of Lystrosaurus we found in the dataset
println("Found ", nrow(lystro_raw), " Lystrosaurus fossils!")
# Found 124 Lystrosaurus fossils!

clean = dropmissing(lystro_raw, [:lng, :lat, :min_ma, :max_ma])
# 124x29 DataFrame; 20 columns and 118 rows omitted

clean = transform(
    clean,
    [:min_ma, :max_ma] => ((lo, hi) -> (lo .+ hi) ./ 2) => :mid_ma
)
# 124x30 DataFrame; 21 columns and 118 rows omitted (This adds :mid_ma to the list of columns)

# Print out how many fossil occurrences of Lystrosaurus have been cleaned up
println("Records with valid coordinates: ", nrow(clean))
println("Records dropped: ", nrow(lystro_raw) - nrow(clean))
# Records with valid coordinates: 124
# Records dropped: 0

# Create a map of Lystrosaurus fossil occurrences
fig = Figure(size = (1000, 600))
ga = GeoAxis(
    fig[1,1],
    dest = "proj=natearth2",
    title = "Lystrosaurus fossil occurrences"
)

poly!(ga, GeoMakie.land(); color = :whitesmoke, strokecolor = :gray55, strokewidth = 0.5)

sc = scatter!(ga, clean.lng, clean.lat; color = clean.mid_ma, colormap = :thermal, markersize = 7.5)

Colorbar(fig[1, 2], sc; label = "Age (Ma)", flipaxis = false)

display(fig)

save("lystrosaurus_pangea_puzzle_occurrences_072326.png", fig)

# This creates a map of Lystrosaurus occurrences, with the age of the fossil (in mya) represented by each point.
# In this case, Lystrosaurus fossils were mainly found in South Africa, but also through Asia and Antarctica.

# Reconstruct the map using a different geological period (e.g., Triassic)

animal = "Lystrosaurus"
period = "Early Triassic"

lystro_triassic_raw = pbdb_occurrences(
    base_name = animal,
    interval = period,
    show = ["coords", "classext"],
    vocab = "pbdb",
    extids = true
)
# 117x29 DataFrame for Triassic Lystrosaurus fossils; 22 columns and 111 rows omitted

println("Found ", nrow(lystro_triassic_raw), " Early Triassic Lystrosaurus fossils!")
# Found 117 Early Triassic Lystrosaurus fossils!

clean = dropmissing(lystro_triassic_raw, [:lng, :lat, :min_ma, :max_ma])

clean = transform(
    clean,
    [:min_ma, :max_ma] => ((lo, hi) -> (lo .+ hi) ./ 2) => :mid_ma
)

println("Records with valid coordinates: ", nrow(clean))
println("Records dropped: ", nrow(lystro_triassic_raw) - nrow(clean))
# Records with valid coordinates: 117
# Records dropped: 0

fig = Figure(size = (1000, 600))
ga = GeoAxis(
    fig[1,1],
    dest = "proj=natearth2",
    title = "Triassic Lystrosaurus fossil occurrences"
)

poly!(ga, GeoMakie.land(); color = :whitesmoke, strokecolor = :gray55, strokewidth = 0.5)

sc = scatter!(ga, clean.lng, clean.lat; color = clean.mid_ma, colormap = :acton, markersize = 7.5)

Colorbar(fig[1, 2], sc; label = "Age (Ma)", flipaxis = false)

display(fig)

save("early_triassic_lystrosaurus_pangea_puzzle_occurrences_072326.png", fig)

# ------------------------------------------

# Part 2: Mesosaurus, and Glossopteris

animal = "Mesosaurus"

mesosaurus_raw = pbdb_occurrences(
    base_name = animal,
    show = ["coords", "classext"],
    vocab = "pbdb",
    extids = true
)
# 30x29 DataFrame for Mesosaurus fossils; 20 columns and 24 rows omitted

println("Found ", nrow(mesosaurus_raw), " Mesosaurus fossils!")
# Found 30 Mesosaurus fossils!

clean = dropmissing(mesosaurus_raw, [:lng, :lat, :min_ma, :max_ma])

clean = transform(
    clean,
    [:min_ma, :max_ma] => ((lo, hi) -> (lo .+ hi) ./ 2) => :mid_ma
)

println("Records with valid coordinates: ", nrow(clean))
println("Records dropped: ", nrow(mesosaurus_raw) - nrow(clean))
# Records with valid coordinates: 30
# Records dropped: 0

fig = Figure(size = (1000, 600))
ga = GeoAxis(
    fig[1,1],
    dest = "proj=natearth2",
    title = "Mesosaurus fossil occurrences"
)

poly!(ga, GeoMakie.land(); color = :whitesmoke, strokecolor = :gray60, strokewidth = 0.75)

sc = scatter!(ga, clean.lng, clean.lat; color = clean.mid_ma, colormap = :lapaz, markersize = 6.5)

Colorbar(fig[1, 2], sc; label = "Age (Ma)", flipaxis = false)

display(fig)

save("mesosaurus_pangea_puzzle_occurrences_072326.png", fig)

# Mesosaurus fossils were mainly distribued throughout South America and Africa, according to the map.
# Nextly, we aim to look up information about Mesosaurus itself, using the function "pbdb_occurrences."

mesosaurus_occs = pbdb_occurrences(
    base_name = "Mesosaurus",
    show = "full",
    vocab = "pbdb",
    extids = true,
)
# 30x134 DataFrame; 125 columns and 24 rows omitted

typeof(mesosaurus_occs)
# This is a DataFrame.

size(mesosaurus_occs)
# (30, 134)

foreach(println, names(mesosaurus_occs))
# Using "foreach" to print out the names of all the columns in the DataFrame, there are 134 columns, including:
# occurrence_no, flags, accepted_name, early_interval, late_interval, lng, lat, min_ma, max_ma, etc.

describe(mesosaurus_occs)
# 134x7 DataFrame; Using "describe" to summarize the DataFrame, there are 134 variables, with numerous statistics
# for each variable, such as mean, minimum, maximum, and nmissing.

show(describe(mesosaurus_occs), allrows = true)
# Using "show" to pull up all variables (as rows), the variables are the same as those in the initial DataFrame.

combine(groupby(mesosaurus_occs, :accepted_rank), nrow)
# 2x2 DataFrame; "Species" and "genus" are considered to be the two accepted ranks for Mesosaurus fossils; nrow = 29 in total.

# If you want to create a visualization of Mesosaurus species, you can use PaleobiologyDB's PhyloPicMakie package, which provides
# a reliable way to visualize the phylogenetic relationships of different fossil organisms using silhouettes like this:
using CairoMakie, PaleobiologyDB, PaleobiologyDB.Taxonomy, PaleobiologyDB.PhyloPicMakie

mesosaurus_species = Taxonomy.child_taxa("Mesosaurus", "species")
# 3-element Vector {String}, which includes the following species:
# "Mesosaurus braziliensis"
# "Mesosaurus capensis"
# "Mesosaurus tenuidens"

phylopic_thumbnail_grid(mesosaurus_species)
# This generates a group of silhouettes for Mesosaurus species that are captured as a three-element Vector (in String)

# Now... when did Mesosaurus live? We can calculate the range of Mesosaurus fossils using occurrences, as follows.

using PaleobiologyDB, PaleobiologyDB.Taxonomy
using DataFrames, Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

animal = "Mesosaurus"

df_meso_raw = pbdb_occurrences(
    base_name = animal,
    show = "full",
    vocab = "pbdb",
    extids = true
)
# You can ensure that you can name your data using unique names, such as "df_meso_raw" for example.
# 30x134 DataFrame; 125 columns and 24 rows omitted

# Data cleaning
df_meso = dropmissing(df_meso_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])
# 30x134 DataFrame (same as raw DataFrame)

# Drop unresolved or unqualified taxa using the following function:
df_meso = drop_unresolved_taxa(df_meso, :genus)
# There was likely no unresolved taxa in this case, as the DataFrame is the same as before (30x134).

# Compute a midpoint age (Ma)
df_meso.midpoint_age = (df_meso.max_ma .+ df_meso.min_ma) ./ 2
# This forms a 30-element Vector{Float64} for each midpoint age (Ma) of the Mesosaurus fossils.

# Using the split-apply-combine method...
appearance_df_meso = combine(
    groupby(df_meso, :genus, sort = true),
    :midpoint_age => maximum => :first_appearance_datum,
    :midpoint_age => minimum => :last_appearance_datum
)
appearance_df_meso = sort(appearance_df_meso, :first_appearance_datum, rev = true)
# 1x3 DataFrame; This shows the first and last appearance datum for the genus of Mesosaurus fossils.

# To build a chart of this particular range...
using Makie, CairoMakie
fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel = "Midpoint age (Ma)",
    ylabel = "Genus",
    yticks = (1:nrow(appearance_df_meso), appearance_df_meso.genus),
    xreversed = true
)

for i in 1:nrow(appearance_df_meso)
    x1 = appearance_df_meso.first_appearance_datum[i]
    x2 = appearance_df_meso.last_appearance_datum[i]
    lines!(ax, [x1, x2], [i, i])
end

fig

save("mesosaurus_range_through_072326.png", fig)

# Repeat the process for Glossopteris
# Overview of Glossopteris fossils

animal = "Glossopteris"

glossopteris_raw = pbdb_occurrences(
    base_name = animal,
    show = ["coords", "classext"],
    vocab = "pbdb",
    extids = true
)

println("Found ", nrow(glossopteris_raw), " Glossopteris fossils!")

clean = dropmissing(glossopteris_raw, [:lng, :lat, :min_ma, :max_ma])

clean = transform(
    clean,
    [:min_ma, :max_ma] => ((lo, hi) -> (lo .+ hi) ./ 2) => :mid_ma
)

println("Records with valid coordinates: ", nrow(clean))
println("Records dropped: ", nrow(glossopteris_raw) - nrow(clean))

fig = Figure(size = (1000, 600))
ga = GeoAxis(
    fig[1,1],
    dest = "proj=natearth2",
    title = "Glossopteris fossil occurrences"
)

poly!(ga, GeoMakie.land(); color = :whitesmoke, strokecolor = :gray60, strokewidth = 0.75)

sc = scatter!(ga, clean.lng, clean.lat; color = clean.mid_ma, colormap = :thermal, markersize = 6.5)

Colorbar(fig[1, 2], sc; label = "Age (Ma)", flipaxis = false)

display(fig)

save("glossopteris_pangea_puzzle_occurrences_072326.png", fig)

# When did Glossopteris live?
using PaleobiologyDB, PaleobiologyDB.Taxonomy
using DataFrames, Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

animal = "Glossopteris"

df_glossopteris_raw = pbdb_occurrences(
    base_name = animal,
    show = "full",
    vocab = "pbdb",
    extids = true
)

df_glossopteris = dropmissing(df_glossopteris_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_glossopteris = drop_unresolved_taxa(df_glossopteris, :genus)

df_glossopteris.midpoint_age = (df_glossopteris.max_ma .+ df_glossopteris.min_ma) ./ 2

appearance_df_glossopteris = combine(
    groupby(df_glossopteris, :genus, sort = true),
    :midpoint_age => maximum => :first_appearance_datum,
    :midpoint_age => minimum => :last_appearance_datum
)
appearance_df_glossopteris = sort(appearance_df_glossopteris, :first_appearance_datum, rev = true)

using Makie, CairoMakie
fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel = "Midpoint age (Ma) of Glossopteris fossils",
    ylabel = "Genus",
    yticks = (1:nrow(appearance_df_glossopteris), appearance_df_glossopteris.genus),
    xreversed = true
)

for i in 1:nrow(appearance_df_glossopteris)
    x1 = appearance_df_glossopteris.first_appearance_datum[i]
    x2 = appearance_df_glossopteris.last_appearance_datum[i]
    lines!(ax, [x1, x2], [i, i])
end

fig

save("glossopteris_range_through_072326.png", fig)

# ------------------------------------------

# Part 3: Constructing maps of Marsupial fossils in the Neogene time period and other earlier time periods
# Then, compare all maps to one another and analyze distribution patterns

# Example: In the Neogene...

animal = "Marsupialia"
period = "Neogene"

marsupialia_neogene_raw = pbdb_occurrences(
    base_name = animal,
    interval = period,
    show = ["coords", "classext"],
    vocab = "pbdb",
    extids = true
)

println("Found ", nrow(marsupialia_neogene_raw), " Neogene Marsupialia fossils!")

clean_marsu = dropmissing(marsupialia_neogene_raw, [:lng, :lat, :min_ma, :max_ma])
clean_marsu = transform(
    clean_marsu,
    [:min_ma, :max_ma] => ((lo, hi) -> (lo .+ hi) ./ 2) => :mid_ma
)

println("Records with valid coordinates: ", nrow(clean_marsu))
println("Records dropped: ", nrow(marsupialia_neogene_raw) - nrow(clean_marsu))

fig = Figure(size = (1000, 600))
ga = GeoAxis(
    fig[1,1],
    dest = "proj=natearth2",
    title = "Marsupialia fossil occurrences in the Neogene"
)

poly!(ga, GeoMakie.land(); color = :whitesmoke, strokecolor = :gray55, strokewidth = 0.5)

sc = scatter!(ga, clean_marsu.lng, clean_marsu.lat; color = clean_marsu.mid_ma, colormap = :thermal, markersize = 7.5)

Colorbar(fig[1, 2], sc; label = "Age (Ma)", flipaxis = false)

display(fig)

save("marsupialia_neogene_occurrences_072326.png", fig)

# Part 4: Choose your own fossil species and repeat the process! What did you learn about the distribution of your chosen fossil species?
# How does it compare to the distribution of other fossils species, such as Lystrosaurus and Marsupialia?