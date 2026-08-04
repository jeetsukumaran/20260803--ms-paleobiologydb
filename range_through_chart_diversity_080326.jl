using PaleobiologyDB, DataFrames, CairoMakie, GeoMakie

# Part 1: Creating a range-through diversity chart for a set of fossils
# Using Organism #1: Edrioasteroidea (class)
focal_name = "Edrioasteroidea"

taxon_names_df = pbdb_taxa(
    name = "Edrioasteroidea",
    rel = "all_children",
    show = "full",
    vocab = "pbdb",
)
# 199x52 DataFrame

taxon_names_df

filter(row -> row.accepted_rank == "family", taxon_names_df)
# 15x52 DataFrame; Mainly including families

occs_edri = pbdb_occurrences(
    base_name = "Edrioasteroidea",
    show = "full",
    vocab = "pbdb"
)
# 289x135 DataFrame

show(occs_edri, allrows = true, allcols = true)
# This gives a "spatial," almost zero-gravity like textual diagram of all names for rows and columns of the DataFrame.

describe(occs_edri)
# 135x7 DataFrame; Consists of statistics such as mean, median, and max/min

filter(row -> row.accepted_rank == "genus", occs_edri)
# 145x135 DataFrame; 145 "names" are all genera

filter(row -> row.accepted_rank == "species", occs_edri)
# 115x135 DataFrame; 145 "names" are all named species

# Now, in full, with genera...
filter(row -> row.accepted_rank == "genus", occs_edri)
occs_edri_genera = filter(row -> !ismissing(row.genus) && !isempty(row.genus), occs_edri)
# 260x135 DataFrame

# Some more analysis

occs_edri
# 289x135 DataFrame

occs_edri[1,1]
# 2262

occs_edri[1, [:phylum, :class, :order, :genus, :accepted_name]]
# One row, which goes from phylum "Echinodermata" to species name "Sprinkleoglobus lloydi"

x = :genus

occs_edri[1, [x]]
# DataFrameRow

# Row │ genus           
#     │ String31        
#─────┼─────────────────
#   1 │ Sprinkleoglobus

# Evaluating ranges from occurrences
using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

# Automatic data caching
PaleobiologyDB.set_autocaching!(true);

occs_edri_raw = pbdb_occurrences(
    base_name = "Edrioasteroidea",
    show = "full",
    vocab = "pbdb"
)
# 289x135 DataFrame

occs_edri = dropmissing(occs_edri_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])
# 289x135 DataFrame

occs_edri = drop_unresolved_taxa(occs_edri, :genus)
# 260x135 DataFrame

occs_edri.midpoint_age = (occs_edri.max_ma .+ occs_edri.min_ma) ./ 2
# 260-element Vector {Float64}

occs_edri = sort(occs_edri, [:genus, :midpoint_age])
# 260x136 DataFrame; New row for Midpoint age (Ma)

# Using the Split-apply-combine
appearance_occs_edri = combine(
    groupby(occs_edri, :genus, sort = true),
    :midpoint_age => maximum => :first_appearance_datum,
    :midpoint_age => minimum => :last_appearance_datum
)
# 60x3 DataFrame; First row consists of genera names; next two raws are first and last appearances respectively

appearance_occs_edri = sort(appearance_occs_edri, :first_appearance_datum, rev=true)
# 60x3 DataFrame

# To manually create a DataFrame for plotting...

range_through_chart_edri = DataFrame(
    genus = appearance_occs_edri.genus,
    first_appearance_datum = appearance_occs_edri.first_appearance_datum,
    last_appearance_datum = appearance_occs_edri.last_appearance_datum
)

using Makie, CairoMakie, ColorSchemes

fig = Figure(size = (800, 700))
ax = Axis(
    fig[1,1],
    xlabel = "Midpoint age (Ma) of Edrioasteroidea genera",
    ylabel = "Genera",
    yticks = (1:nrow(appearance_occs_edri), appearance_occs_edri.genus),
    yticklabelsize = 7,
    title = "Ranges and diversity of Edrioasteroidea genera",
    xreversed = true
)

for i in 1:nrow(appearance_occs_edri)
    x1 = appearance_occs_edri.first_appearance_datum[i]
    x2 = appearance_occs_edri.last_appearance_datum[i]
    lines!(ax, [x1, x2], [i, i], colorscale = :auto)
end

fig

# Part 2: Creating a diversity curve
# In this scenario, you'll want to create a boundary results table with a diversity loop using interval boundaries

using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

# This time, you will set up auto-caching
PaleobiologyDB.set_autocaching!(true);

# Again, using Edrioasteroidea as the example fossil group of interest
df_edri_diversity_raw = pbdb_occurrences(
    base_name = "Edrioasteroidea",
    show = "full",
    vocab = "pbdb"
)

# Data cleaning
df_edri_diversity = dropmissing(df_edri_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_edri_diversity = Taxonomy.drop_unqualified_taxa(df_edri_diversity, "genus")

df_edri_diversity.midpoint_age = (df_edri_diversity.max_ma .+ df_edri_diversity.min_ma) ./ 2

interval_df_1 = combine(
    groupby(df_edri_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_1.interval_mid_ma = (
    interval_df_1.interval_max_ma .+ interval_df_1.interval_min_ma
) ./2

interval_df_1 = sort(interval_df_1, :interval_mid_ma, rev=true)

genus_ranges_of_edri = combine(
    groupby(df_edri_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

# Creating a diversity loop through time
diversity = Int[]

for interval_row in eachrow(interval_df_1)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_edri)
    )
    push!(diversity, n)
end

interval_df_1.diversity = diversity

# For the nineteen other fossils...

using PaleobiologyDB, PaleobiologyDB.Taxonomy, DataFrames
using Makie, CairoMakie, GeoMakie, ColorSchemes

PaleobiologyDB.set_autocaching!(true);

## (1) Articulata

df_arti_diversity_raw = pbdb_occurrences(
    base_name = "Articulata",
    show = "full",
    vocab = "pbdb"
)

df_arti_diversity = dropmissing(df_arti_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_arti_diversity = Taxonomy.drop_unqualified_taxa(df_arti_diversity, "genus")

df_arti_diversity.midpoint_age = (df_arti_diversity.max_ma .+ df_arti_diversity.min_ma) ./ 2

interval_df_2 = combine(
    groupby(df_arti_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_2.interval_mid_ma = (
    interval_df_2.interval_max_ma .+ interval_df_2.interval_min_ma
) ./2

interval_df_2 = sort(interval_df_2, :interval_mid_ma, rev=true)

genus_ranges_of_arti = combine(
    groupby(df_arti_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_2)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_arti)
    )
    push!(diversity, n)
end

interval_df_2.diversity = diversity

## (2) Ammonoidea

df_ammo_diversity_raw = pbdb_occurrences(
    base_name = "Ammonoidea",
    show = "full",
    vocab = "pbdb"
)

df_ammo_diversity = dropmissing(df_ammo_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_ammo_diversity = Taxonomy.drop_unqualified_taxa(df_ammo_diversity, "genus")

df_ammo_diversity.midpoint_age = (df_ammo_diversity.max_ma .+ df_ammo_diversity.min_ma) ./ 2

interval_df_3 = combine(
    groupby(df_ammo_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_3.interval_mid_ma = (
    interval_df_3.interval_max_ma .+ interval_df_3.interval_min_ma
) ./2

interval_df_3 = sort(interval_df_3, :interval_mid_ma, rev=true)

genus_ranges_of_ammo = combine(
    groupby(df_ammo_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_3)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_ammo)
    )
    push!(diversity, n)
end

interval_df_3.diversity = diversity

## (3) Rugosa

df_rugosa_diversity_raw = pbdb_occurrences(
    base_name = "Rugosa",
    show = "full",
    vocab = "pbdb"
)

df_rugosa_diversity = dropmissing(df_rugosa_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_rugosa_diversity = Taxonomy.drop_unqualified_taxa(df_rugosa_diversity, "genus")

df_rugosa_diversity.midpoint_age = (df_rugosa_diversity.max_ma .+ df_rugosa_diversity.min_ma) ./ 2

interval_df_4 = combine(
    groupby(df_rugosa_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_4.interval_mid_ma = (
    interval_df_4.interval_max_ma .+ interval_df_4.interval_min_ma
) ./2

interval_df_4 = sort(interval_df_4, :interval_mid_ma, rev=true)

genus_ranges_of_rugosa = combine(
    groupby(df_rugosa_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_4)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_rugosa)
    )
    push!(diversity, n)
end

interval_df_4.diversity = diversity

## (4) Tabulata

df_tabu_diversity_raw = pbdb_occurrences(
    base_name = "Tabulata",
    show = "full",
    vocab = "pbdb"
)

df_tabu_diversity = dropmissing(df_tabu_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_tabu_diversity = Taxonomy.drop_unqualified_taxa(df_tabu_diversity, "genus")

df_tabu_diversity.midpoint_age = (df_tabu_diversity.max_ma .+ df_tabu_diversity.min_ma) ./ 2

interval_df_5 = combine(
    groupby(df_tabu_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_5.interval_mid_ma = (
    interval_df_5.interval_max_ma .+ interval_df_5.interval_min_ma
) ./2

interval_df_5 = sort(interval_df_5, :interval_mid_ma, rev=true)

genus_ranges_of_tabu = combine(
    groupby(df_tabu_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_5)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_tabu)
    )
    push!(diversity, n)
end

interval_df_5.diversity = diversity

## (5) Blastoidea

df_blast_diversity_raw = pbdb_occurrences(
    base_name = "Blastoidea",
    show = "full",
    vocab = "pbdb"
)

df_blast_diversity = dropmissing(df_blast_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_blast_diversity = Taxonomy.drop_unqualified_taxa(df_blast_diversity, "genus")

df_blast_diversity.midpoint_age = (df_blast_diversity.max_ma .+ df_blast_diversity.min_ma) ./ 2

interval_df_6 = combine(
    groupby(df_blast_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_6.interval_mid_ma = (
    interval_df_6.interval_max_ma .+ interval_df_6.interval_min_ma
) ./2

interval_df_6 = sort(interval_df_6, :interval_mid_ma, rev=true)

genus_ranges_of_blast = combine(
    groupby(df_blast_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_6)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_blast)
    )
    push!(diversity, n)
end

interval_df_6.diversity = diversity

## (6) Rostroconchia

df_rost_diversity_raw = pbdb_occurrences(
    base_name = "Rostroconchia",
    show = "full",
    vocab = "pbdb"
)

df_rost_diversity = dropmissing(df_rost_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_rost_diversity = Taxonomy.drop_unqualified_taxa(df_rost_diversity, "genus")

df_rost_diversity.midpoint_age = (df_rost_diversity.max_ma .+ df_rost_diversity.min_ma) ./ 2

interval_df_7 = combine(
    groupby(df_rost_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_7.interval_mid_ma = (
    interval_df_7.interval_max_ma .+ interval_df_7.interval_min_ma
) ./2

interval_df_7 = sort(interval_df_7, :interval_mid_ma, rev=true)

genus_ranges_of_rost = combine(
    groupby(df_rost_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_7)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_rost)
    )
    push!(diversity, n)
end

interval_df_7.diversity = diversity

## (7) Gastropoda

df_gast_diversity_raw = pbdb_occurrences(
    base_name = "Gastropoda",
    show = "full",
    vocab = "pbdb"
)

df_gast_diversity = dropmissing(df_gast_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_gast_diversity = Taxonomy.drop_unqualified_taxa(df_gast_diversity, "genus")

df_gast_diversity.midpoint_age = (df_gast_diversity.max_ma .+ df_gast_diversity.min_ma) ./ 2

interval_df_8 = combine(
    groupby(df_gast_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_8.interval_mid_ma = (
    interval_df_8.interval_max_ma .+ interval_df_8.interval_min_ma
) ./2

interval_df_8 = sort(interval_df_8, :interval_mid_ma, rev=true)

genus_ranges_of_gast = combine(
    groupby(df_gast_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_8)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_gast)
    )
    push!(diversity, n)
end

interval_df_8.diversity = diversity

## (8) Stenolaemata

df_sten_diversity_raw = pbdb_occurrences(
    base_name = "Stenolaemata",
    show = "full",
    vocab = "pbdb"
)

df_sten_diversity = dropmissing(df_sten_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_sten_diversity = Taxonomy.drop_unqualified_taxa(df_sten_diversity, "genus")

df_sten_diversity.midpoint_age = (df_sten_diversity.max_ma .+ df_sten_diversity.min_ma) ./ 2

interval_df_9 = combine(
    groupby(df_sten_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_9.interval_mid_ma = (
    interval_df_9.interval_max_ma .+ interval_df_9.interval_min_ma
) ./2

interval_df_9 = sort(interval_df_9, :interval_mid_ma, rev=true)

genus_ranges_of_sten = combine(
    groupby(df_sten_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_9)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_sten)
    )
    push!(diversity, n)
end

interval_df_9.diversity = diversity

## (9) Trilobita

df_tri_diversity_raw = pbdb_occurrences(
    base_name = "Trilobita",
    show = "full",
    vocab = "pbdb"
)

df_tri_diversity = dropmissing(df_tri_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_tri_diversity = Taxonomy.drop_unqualified_taxa(df_tri_diversity, "genus")

df_tri_diversity.midpoint_age = (df_tri_diversity.max_ma .+ df_tri_diversity.min_ma) ./ 2

interval_df_10 = combine(
    groupby(df_tri_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_10.interval_mid_ma = (
    interval_df_10.interval_max_ma .+ interval_df_10.interval_min_ma
) ./2

interval_df_10 = sort(interval_df_10, :interval_mid_ma, rev=true)

genus_ranges_of_tri = combine(
    groupby(df_tri_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_10)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_tri)
    )
    push!(diversity, n)
end

interval_df_10.diversity = diversity

## (10) Scleractinia

df_sclera_diversity_raw = pbdb_occurrences(
    base_name = "Scleractinia",
    show = "full",
    vocab = "pbdb"
)

df_sclera_diversity = dropmissing(df_sclera_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_sclera_diversity = Taxonomy.drop_unqualified_taxa(df_sclera_diversity, "genus")

df_sclera_diversity.midpoint_age = (df_sclera_diversity.max_ma .+ df_sclera_diversity.min_ma) ./ 2

interval_df_11 = combine(
    groupby(df_sclera_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_11.interval_mid_ma = (
    interval_df_11.interval_max_ma .+ interval_df_11.interval_min_ma
) ./2

interval_df_11 = sort(interval_df_11, :interval_mid_ma, rev=true)

genus_ranges_of_sclera = combine(
    groupby(df_sclera_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_11)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_sclera)
    )
    push!(diversity, n)
end

interval_df_11.diversity = diversity

## (11) Crinoidea

df_crino_diversity_raw = pbdb_occurrences(
    base_name = "Crinoidea",
    show = "full",
    vocab = "pbdb"
)

df_crino_diversity = dropmissing(df_crino_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_crino_diversity = Taxonomy.drop_unqualified_taxa(df_crino_diversity, "genus")

df_crino_diversity.midpoint_age = (df_crino_diversity.max_ma .+ df_crino_diversity.min_ma) ./ 2

interval_df_12 = combine(
    groupby(df_crino_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_12.interval_mid_ma = (
    interval_df_12.interval_max_ma .+ interval_df_12.interval_min_ma
) ./2

interval_df_12 = sort(interval_df_12, :interval_mid_ma, rev=true)

genus_ranges_of_crino = combine(
    groupby(df_crino_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_12)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_crino)
    )
    push!(diversity, n)
end

interval_df_12.diversity = diversity

## (12) Echinoidea

df_echin_diversity_raw = pbdb_occurrences(
    base_name = "Echinoidea",
    show = "full",
    vocab = "pbdb"
)

df_echin_diversity = dropmissing(df_echin_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_echin_diversity = Taxonomy.drop_unqualified_taxa(df_echin_diversity, "genus")

df_echin_diversity.midpoint_age = (df_echin_diversity.max_ma .+ df_echin_diversity.min_ma) ./ 2

interval_df_13 = combine(
    groupby(df_echin_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_13.interval_mid_ma = (
    interval_df_13.interval_max_ma .+ interval_df_13.interval_min_ma
) ./2

interval_df_13 = sort(interval_df_13, :interval_mid_ma, rev=true)

genus_ranges_of_echin = combine(
    groupby(df_echin_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_13)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_echin)
    )
    push!(diversity, n)
end

interval_df_13.diversity = diversity

## (13) Bivalvia

df_biva_diversity_raw = pbdb_occurrences(
    base_name = "Bivalvia",
    show = "full",
    vocab = "pbdb"
)

df_biva_diversity = dropmissing(df_biva_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_biva_diversity = Taxonomy.drop_unqualified_taxa(df_biva_diversity, "genus")

df_biva_diversity.midpoint_age = (df_biva_diversity.max_ma .+ df_biva_diversity.min_ma) ./ 2

interval_df_14 = combine(
    groupby(df_biva_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_14.interval_mid_ma = (
    interval_df_14.interval_max_ma .+ interval_df_14.interval_min_ma
) ./2

interval_df_14 = sort(interval_df_14, :interval_mid_ma, rev=true)

genus_ranges_of_biva = combine(
    groupby(df_biva_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_14)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_biva)
    )
    push!(diversity, n)
end

interval_df_14.diversity = diversity

## (14) Nautiloidea

df_naut_diversity_raw = pbdb_occurrences(
    base_name = "Nautiloidea",
    show = "full",
    vocab = "pbdb"
)

df_naut_diversity = dropmissing(df_naut_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_naut_diversity = Taxonomy.drop_unqualified_taxa(df_naut_diversity, "genus")

df_naut_diversity.midpoint_age = (df_naut_diversity.max_ma .+ df_naut_diversity.min_ma) ./ 2

interval_df_15 = combine(
    groupby(df_naut_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_15.interval_mid_ma = (
    interval_df_15.interval_max_ma .+ interval_df_15.interval_min_ma
) ./2

interval_df_15 = sort(interval_df_15, :interval_mid_ma, rev=true)

genus_ranges_of_naut = combine(
    groupby(df_naut_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_15)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_naut)
    )
    push!(diversity, n)
end

interval_df_15.diversity = diversity

## (15) Stromatoporata (Stromatoporoidea)

df_stro_diversity_raw = pbdb_occurrences(
    base_name = "Stromatoporoidea",
    show = "full",
    vocab = "pbdb"
)

df_stro_diversity = dropmissing(df_stro_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_stro_diversity = Taxonomy.drop_unqualified_taxa(df_stro_diversity, "genus")

df_stro_diversity.midpoint_age = (df_stro_diversity.max_ma .+ df_stro_diversity.min_ma) ./ 2

interval_df_16 = combine(
    groupby(df_stro_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_16.interval_mid_ma = (
    interval_df_16.interval_max_ma .+ interval_df_16.interval_min_ma
) ./2

interval_df_16 = sort(interval_df_16, :interval_mid_ma, rev=true)

genus_ranges_of_stro = combine(
    groupby(df_stro_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_16)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_stro)
    )
    push!(diversity, n)
end

interval_df_16.diversity = diversity

## (16) Eurypterida

df_eury_diversity_raw = pbdb_occurrences(
    base_name = "Eurypterida",
    show = "full",
    vocab = "pbdb"
)

df_eury_diversity = dropmissing(df_eury_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_eury_diversity = Taxonomy.drop_unqualified_taxa(df_eury_diversity, "genus")

df_eury_diversity.midpoint_age = (df_eury_diversity.max_ma .+ df_eury_diversity.min_ma) ./ 2

interval_df_17 = combine(
    groupby(df_eury_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_17.interval_mid_ma = (
    interval_df_17.interval_max_ma .+ interval_df_17.interval_min_ma
) ./2

interval_df_17 = sort(interval_df_17, :interval_mid_ma, rev=true)

genus_ranges_of_eury = combine(
    groupby(df_eury_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_17)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_eury)
    )
    push!(diversity, n)
end

interval_df_17.diversity = diversity

## (17) Graptolithina

df_grapt_diversity_raw = pbdb_occurrences(
    base_name = "Graptolithina",
    show = "full",
    vocab = "pbdb"
)

df_grapt_diversity = dropmissing(df_grapt_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_grapt_diversity = Taxonomy.drop_unqualified_taxa(df_grapt_diversity, "genus")

df_grapt_diversity.midpoint_age = (df_grapt_diversity.max_ma .+ df_grapt_diversity.min_ma) ./ 2

interval_df_18 = combine(
    groupby(df_grapt_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_18.interval_mid_ma = (
    interval_df_18.interval_max_ma .+ interval_df_18.interval_min_ma
) ./2

interval_df_18 = sort(interval_df_18, :interval_mid_ma, rev=true)

genus_ranges_of_grapt = combine(
    groupby(df_grapt_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_18)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_grapt)
    )
    push!(diversity, n)
end

interval_df_18.diversity = diversity

## (18) Conodonta

df_cono_diversity_raw = pbdb_occurrences(
    base_name = "Conodonta",
    show = "full",
    vocab = "pbdb"
)

df_cono_diversity = dropmissing(df_cono_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_cono_diversity = Taxonomy.drop_unqualified_taxa(df_cono_diversity, "genus")

df_cono_diversity.midpoint_age = (df_cono_diversity.max_ma .+ df_cono_diversity.min_ma) ./ 2

interval_df_19 = combine(
    groupby(df_cono_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_19.interval_mid_ma = (
    interval_df_19.interval_max_ma .+ interval_df_19.interval_min_ma
) ./2

interval_df_19 = sort(interval_df_19, :interval_mid_ma, rev=true)

genus_ranges_of_cono = combine(
    groupby(df_cono_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_19)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_cono)
    )
    push!(diversity, n)
end

interval_df_19.diversity = diversity

## (19) Insecta

df_insect_diversity_raw = pbdb_occurrences(
    base_name = "Insecta",
    show = "full",
    vocab = "pbdb"
)

df_insect_diversity = dropmissing(df_insect_diversity_raw, [
    :genus,
    :early_interval,
    :max_ma,
    :min_ma
])

df_insect_diversity = Taxonomy.drop_unqualified_taxa(df_insect_diversity, "genus")

df_insect_diversity.midpoint_age = (df_insect_diversity.max_ma .+ df_insect_diversity.min_ma) ./ 2

interval_df_20 = combine(
    groupby(df_insect_diversity, :early_interval),
    :max_ma => maximum => :interval_max_ma,
    :min_ma => minimum => :interval_min_ma,
)

interval_df_20.interval_mid_ma = (
    interval_df_20.interval_max_ma .+ interval_df_20.interval_min_ma
) ./2

interval_df_20 = sort(interval_df_20, :interval_mid_ma, rev=true)

genus_ranges_of_insect = combine(
    groupby(df_insect_diversity, :genus),
    :midpoint_age => maximum => :fad_ma,
    :midpoint_age => minimum => :lad_ma,
)

diversity = Int[]

for interval_row in eachrow(interval_df_20)
    n = count(
        g -> g.fad_ma >= interval_row.interval_mid_ma >= g.lad_ma,
        eachrow(genus_ranges_of_insect)
    )
    push!(diversity, n)
end

interval_df_20.diversity = diversity

# Create DataFrames consisting of both interval and diversity data for Edrioasteroidea and all other fossils; As such:

range_through_1 = DataFrame(
    edrioasteroidea_diversity = interval_df_1.diversity,
    interval_mid_ma_edri = interval_df_1.interval_mid_ma
)

range_through_2 = DataFrame(
    articulata_diversity = interval_df_2.diversity,
    interval_mid_ma_arti = interval_df_2.interval_mid_ma
)

range_through_3 = DataFrame(
    ammonoidea_diversity = interval_df_3.diversity,
    interval_mid_ma_ammo = interval_df_3.interval_mid_ma
)

range_through_4 = DataFrame(
    rugosa_diversity = interval_df_4.diversity,
    interval_mid_ma_rugosa = interval_df_4.interval_mid_ma
)

range_through_5 = DataFrame(
    tabulata_diversity = interval_df_5.diversity,
    interval_mid_ma_tabu = interval_df_5.interval_mid_ma
)

range_through_6 = DataFrame(
    blastoidea_diversity = interval_df_6.diversity,
    interval_mid_ma_blast = interval_df_6.interval_mid_ma
)

range_through_7 = DataFrame(
    rostroconchia_diversity = interval_df_7.diversity,
    interval_mid_ma_rost = interval_df_7.interval_mid_ma
)

range_through_8 = DataFrame(
    gastropoda_diversity = interval_df_8.diversity,
    interval_mid_ma_gast = interval_df_8.interval_mid_ma
)

range_through_9 = DataFrame(
    stenolaemata_diversity = interval_df_9.diversity,
    interval_mid_ma_sten = interval_df_9.interval_mid_ma
)

range_through_10 = DataFrame(
    trilobita_diversity = interval_df_10.diversity,
    interval_mid_ma_tri = interval_df_10.interval_mid_ma
)

range_through_11 = DataFrame(
    scleractinia_diversity = interval_df_11.diversity,
    interval_mid_ma_sclera = interval_df_11.interval_mid_ma
)

range_through_12 = DataFrame(
    crinoidea_diversity = interval_df_12.diversity,
    interval_mid_ma_crino = interval_df_12.interval_mid_ma
)

range_through_13 = DataFrame(
    echinoidea_diversity = interval_df_13.diversity,
    interval_mid_ma_echin = interval_df_13.interval_mid_ma
)

range_through_14 = DataFrame(
    bivalvia_diversity = interval_df_14.diversity,
    interval_mid_ma_biva = interval_df_14.interval_mid_ma
)

range_through_15 = DataFrame(
    nautiloidea_diversity = interval_df_15.diversity,
    interval_mid_ma_naut = interval_df_15.interval_mid_ma
)

range_through_16 = DataFrame(
    stromatoporata_diversity = interval_df_16.diversity,
    interval_mid_ma_stro = interval_df_16.interval_mid_ma
)

range_through_17 = DataFrame(
    eurypterida_diversity = interval_df_17.diversity,
    interval_mid_ma_eury = interval_df_17.interval_mid_ma
)

range_through_18 = DataFrame(
    graptolithina_diversity = interval_df_18.diversity,
    interval_mid_ma_grapt = interval_df_18.interval_mid_ma
)

range_through_19 = DataFrame(
    conodonta_diversity = interval_df_19.diversity,
    interval_mid_ma_cono = interval_df_19.interval_mid_ma
)

range_through_20 = DataFrame(
    insecta_diversity = interval_df_20.diversity,
    interval_mid_ma_insect = interval_df_20.interval_mid_ma
)

# Creating diversity curves with our DataFrames

# Curve 1)

using Makie, CairoMakie, ColorSchemes
fig = Figure(size = (1200, 900), fontsize = 11)
ax1 = Axis(
    fig[1, 1],
    xlabel = "Age (Ma)",
    ylabel = "Number of genera",
    title = "Genus diversity of several fossils through time",
    xticks = (interval_df_1.interval_mid_ma, interval_df_1.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax2 = Axis(
    fig[1,1],
    xticks = (interval_df_2.interval_mid_ma, interval_df_2.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax3 = Axis(
    fig[1,1],
    xticks = (interval_df_3.interval_mid_ma, interval_df_3.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax4 = Axis(
    fig[1,1],
    xticks = (interval_df_4.interval_mid_ma, interval_df_4.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax5 = Axis(
    fig[1,1],
    xticks = (interval_df_5.interval_mid_ma, interval_df_5.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
hidedecorations!(ax2)
hidedecorations!(ax3)
hidedecorations!(ax4)
hidedecorations!(ax5)

lines!(ax1, interval_df_1.interval_mid_ma, interval_df_1.diversity, linewidth = 2, color = :red4)
scatter!(ax1, interval_df_1.interval_mid_ma, interval_df_1.diversity, markersize = 6, color = :red4, label = "Edrioasteroidea")

lines!(ax2, interval_df_2.interval_mid_ma, interval_df_2.diversity, linewidth = 2, color = :goldenrod)
scatter!(ax2, interval_df_2.interval_mid_ma, interval_df_2.diversity, markersize = 6, color = :goldenrod, label = "Articulata")

lines!(ax3, interval_df_3.interval_mid_ma, interval_df_3.diversity, linewidth = 2, color = :dodgerblue)
scatter!(ax3, interval_df_3.interval_mid_ma, interval_df_3.diversity, markersize = 6, color = :dodgerblue, label = "Ammonoidea")

lines!(ax4, interval_df_4.interval_mid_ma, interval_df_4.diversity, linewidth = 2, color = :palevioletred)
scatter!(ax4, interval_df_4.interval_mid_ma, interval_df_4.diversity, markersize = 6, color = :palevioletred, label = "Rugosa")

lines!(ax5, interval_df_5.interval_mid_ma, interval_df_5.diversity, linewidth = 2, color = :darkorchid3)
scatter!(ax5, interval_df_5.interval_mid_ma, interval_df_5.diversity, markersize = 6, color = :darkorchid3, label = "Tabulata")

axislegend(ax1, position = (0.01, 0.975), framevisible = true, labelsize = 11)
axislegend(ax2, position = (0.01, 0.925), framevisible = true, labelsize = 11)
axislegend(ax3, position = (0.01, 0.875), framevisible = true, labelsize = 11)
axislegend(ax4, position = (0.01, 0.825), framevisible = true, labelsize = 11)
axislegend(ax5, position = (0.01, 0.775), framevisible = true, labelsize = 11)

fig

# Curve 2)

fig = Figure(size = (1200, 900), fontsize = 11)
ax6 = Axis(
    fig[1,2],
    xlabel = "Age (Ma)",
    ylabel = "Number of genera",
    title = "Genus diversity of several fossils through time",
    xticks = (interval_df_6.interval_mid_ma, interval_df_6.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax7 = Axis(
    fig[1,2],
    xticks = (interval_df_7.interval_mid_ma, interval_df_7.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax8 = Axis(
    fig[1,2],
    xticks = (interval_df_8.interval_mid_ma, interval_df_8.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax9 = Axis(
    fig[1,2],
    xticks = (interval_df_9.interval_mid_ma, interval_df_9.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax10 = Axis(
    fig[1,2],
    xticks = (interval_df_10.interval_mid_ma, interval_df_10.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
hidedecorations!(ax7)
hidedecorations!(ax8)
hidedecorations!(ax9)
hidedecorations!(ax10)

lines!(ax6, interval_df_6.interval_mid_ma, interval_df_6.diversity, linewidth = 2, color = :sienna)
scatter!(ax6, interval_df_6.interval_mid_ma, interval_df_6.diversity, markersize = 6, color = :sienna, label = "Blastoidea")

lines!(ax7, interval_df_7.interval_mid_ma, interval_df_7.diversity, linewidth = 2, color = :mediumturquoise)
scatter!(ax7, interval_df_7.interval_mid_ma, interval_df_7.diversity, markersize = 6, color = :mediumturquoise, label = "Rostroconchia")

lines!(ax8, interval_df_8.interval_mid_ma, interval_df_8.diversity, linewidth = 2, color = :ivory4)
scatter!(ax8, interval_df_8.interval_mid_ma, interval_df_8.diversity, markersize = 6, color = :ivory4, label = "Gastropoda")

lines!(ax9, interval_df_9.interval_mid_ma, interval_df_9.diversity, linewidth = 2, color = :coral2)
scatter!(ax9, interval_df_9.interval_mid_ma, interval_df_9.diversity, markersize = 6, color = :coral2, label = "Stenolaemata")

lines!(ax10, interval_df_10.interval_mid_ma, interval_df_10.diversity, linewidth = 2, color = :chartreuse2)
scatter!(ax10, interval_df_10.interval_mid_ma, interval_df_10.diversity, markersize = 6, color = :chartreuse2, label = "Trilobita")

axislegend(ax6, position = (0.01, 0.975), framevisible = true, labelsize = 11)
axislegend(ax7, position = (0.01, 0.925), framevisible = true, labelsize = 11)
axislegend(ax8, position = (0.01, 0.875), framevisible = true, labelsize = 11)
axislegend(ax9, position = (0.01, 0.825), framevisible = true, labelsize = 11)
axislegend(ax10, position = (0.01, 0.775), framevisible = true, labelsize = 11)

fig

# Curve 3)

fig = Figure(size = (1200, 900), fontsize = 11)
ax11 = Axis(
    fig[2,1],
    xlabel = "Age (Ma)",
    ylabel = "Number of genera",
    title = "Genus diversity of several fossils through time",
    xticks = (interval_df_11.interval_mid_ma, interval_df_11.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax12 = Axis(
    fig[2,1],
    xticks = (interval_df_12.interval_mid_ma, interval_df_12.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax13 = Axis(
    fig[2,1],
    xticks = (interval_df_13.interval_mid_ma, interval_df_13.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax14 = Axis(
    fig[2,1],
    xticks = (interval_df_14.interval_mid_ma, interval_df_14.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax15 = Axis(
    fig[2,1],
    xticks = (interval_df_15.interval_mid_ma, interval_df_15.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
hidedecorations!(ax12)
hidedecorations!(ax13)
hidedecorations!(ax14)
hidedecorations!(ax15)

lines!(ax11, interval_df_11.interval_mid_ma, interval_df_11.diversity, linewidth = 2, color = :mediumslateblue)
scatter!(ax11, interval_df_11.interval_mid_ma, interval_df_11.diversity, markersize = 6, color = :mediumslateblue, label = "Scleractinia")

lines!(ax12, interval_df_12.interval_mid_ma, interval_df_12.diversity, linewidth = 2, color = :maroon3)
scatter!(ax12, interval_df_12.interval_mid_ma, interval_df_12.diversity, markersize = 6, color = :maroon3, label = "Crinoidea")

lines!(ax13, interval_df_13.interval_mid_ma, interval_df_13.diversity, linewidth = 2, color = :plum4)
scatter!(ax13, interval_df_13.interval_mid_ma, interval_df_13.diversity, markersize = 6, color = :plum4, label = "Echinoidea")

lines!(ax14, interval_df_14.interval_mid_ma, interval_df_14.diversity, linewidth = 2, color = :darkorange)
scatter!(ax14, interval_df_14.interval_mid_ma, interval_df_14.diversity, markersize = 6, color = :darkorange, label = "Bivalvia")

lines!(ax15, interval_df_15.interval_mid_ma, interval_df_15.diversity, linewidth = 2, color = :blue)
scatter!(ax15, interval_df_15.interval_mid_ma, interval_df_15.diversity, markersize = 6, color = :blue, label = "Nautiloidea")

axislegend(ax11, position = (0.01, 0.975), framevisible = true, labelsize = 11)
axislegend(ax12, position = (0.01, 0.925), framevisible = true, labelsize = 11)
axislegend(ax13, position = (0.01, 0.875), framevisible = true, labelsize = 11)
axislegend(ax14, position = (0.01, 0.825), framevisible = true, labelsize = 11)
axislegend(ax15, position = (0.01, 0.775), framevisible = true, labelsize = 11)

fig

# Curve 4)

fig = Figure(size = (1200, 900), fontsize = 11)
ax16 = Axis(
    fig[2,2],
    xlabel = "Age (Ma)",
    ylabel = "Number of genera",
    title = "Genus diversity of several fossils through time",
    xticks = (interval_df_16.interval_mid_ma, interval_df_16.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax17 = Axis(
    fig[2,2],
    xticks = (interval_df_17.interval_mid_ma, interval_df_17.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax18 = Axis(
    fig[2,2],
    xticks = (interval_df_18.interval_mid_ma, interval_df_18.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax19 = Axis(
    fig[2,2],
    xticks = (interval_df_19.interval_mid_ma, interval_df_19.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
ax20 = Axis(
    fig[2,2],
    xticks = (interval_df_20.interval_mid_ma, interval_df_20.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
hidedecorations!(ax17)
hidedecorations!(ax18)
hidedecorations!(ax19)
hidedecorations!(ax20)

lines!(ax16, interval_df_16.interval_mid_ma, interval_df_16.diversity, linewidth = 2, color = :chocolate2)
scatter!(ax16, interval_df_16.interval_mid_ma, interval_df_16.diversity, markersize = 6, color = :chocolate2, label = "Stromatoporata")

lines!(ax17, interval_df_17.interval_mid_ma, interval_df_17.diversity, linewidth = 2, color = :gray21)
scatter!(ax17, interval_df_17.interval_mid_ma, interval_df_17.diversity, markersize = 6, color = :gray21, label = "Eurypterida")

lines!(ax18, interval_df_18.interval_mid_ma, interval_df_18.diversity, linewidth = 2, color = :darkkhaki)
scatter!(ax18, interval_df_18.interval_mid_ma, interval_df_18.diversity, markersize = 6, color = :darkkhaki, label = "Graptolithina")

lines!(ax19, interval_df_19.interval_mid_ma, interval_df_19.diversity, linewidth = 2, color = :red)
scatter!(ax19, interval_df_19.interval_mid_ma, interval_df_19.diversity, markersize = 6, color = :red, label = "Conodonta")

lines!(ax20, interval_df_20.interval_mid_ma, interval_df_20.diversity, linewidth = 2, color = :green4)
scatter!(ax20, interval_df_20.interval_mid_ma, interval_df_20.diversity, markersize = 6, color = :green4, label = "Insecta")

axislegend(ax16, position = (0.01, 0.975), framevisible = true, labelsize = 11)
axislegend(ax17, position = (0.01, 0.925), framevisible = true, labelsize = 11)
axislegend(ax18, position = (0.01, 0.875), framevisible = true, labelsize = 11)
axislegend(ax19, position = (0.01, 0.825), framevisible = true, labelsize = 11)
axislegend(ax20, position = (0.01, 0.775), framevisible = true, labelsize = 11)

fig

# Currently works for one axis, as demonstrated below in another example; However, it needs more work to fit
# multiple axes (shown above)

# To make this complex code even simpler...

using Makie, CairoMakie, ColorSchemes

fig = Figure(fontsize = 11)
ax = Axis(
    fig[1, 1],
    xlabel = "Age (Ma)",
    ylabel = "Number of genera",
    title = "Genus diversity of Edrioasteroidea fossils through time",
    xticks = (interval_df_1.interval_mid_ma, interval_df_1.early_interval),
    xticklabelrotation = π/4,
    xticklabelsize = 6.5,
    xreversed = true,
)
lines!(ax, interval_df_1.interval_mid_ma, interval_df_1.diversity, linewidth = 2, color = :green)
scatter!(ax, interval_df_1.interval_mid_ma, interval_df_1.diversity, markersize = 8, color = :green, label = "Edrioasteroidea")

axislegend(position = :lt)

fig

# Answer additional questions and analyze data; determine which extinction event was the biggest 