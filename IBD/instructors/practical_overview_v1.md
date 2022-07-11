# Practical Overview
- *Duration*: 1.5 hrs


## Key Concepts
- IBS is affected by PLAF and can produce biased estimates of relatedness
- IBD theoretically follows isolation by distance; when it doesn't, interesting things are happening
- IBD captures recent relatedness, useful for determining transmission events/transmission intensity, importation, corridors of gene flow
  + Therefore, IBS/IBD can measure connectivity on spatial/temporal scales relevant to control

## Part 00: Simulate Data (for instructors to do prior to practical)
- Use `polySimIBD`[https://github.com/nickbrazeau/polySimIBD] to simulate 5 demes of varying transmission intensity with one deme as "island"
    + Layer in IBS manually into the VCF w/ fixed mutation function (or just brute it)?

## Part 0: Set up environment, Read in Data (10 minutes)
- Install multiple packages: source Verity MLE from R script (versus downloading MIPAnalyzer and dealing w/ biallelic MIP class)
  + vcfR?
  + sf?
  + tidygraph?

## Part 1: Identity by State (20 minutes)
- Visualize genotypes
  + Can participants visualize patterns of IBS vs IBD (_not really_!)
- Participants write IBS function
  + Demonstrate that different allele frequency distributions lead to very different mean IBS
- Make histogram of relatedness by IBS --> same or different tails as IBD (_different_)

## Part 2: Identity by Descent (20 minutes)
- Run IBD Calculations
  + Run Verity MLE inbreeding function
  + Run hmmIBD (use Rcpp wrapper from OJ)
- Make histogram of relatedness by IBD --> same or different tails as IBS (_different_)
  + Sites w/ diverse PLAF informative for IBD -> find and show?

## Part 3: Connectivity, Transmission Chains, Importation (30 minutes)
- Isolation by distance: maps
- Networks: between and within deme connections (rule of transitivity)
- Contrast w/ measures of population differentiation (vs individualistic pairwise IBD)
- Look at COI vs IBD segment length --> transmission intensity
  + Can hmmIBD give us segment lengths, no right? can do proportion of highly related pairs instead

## Part 4: Questions & Advanced Section
- Polyclonal infections
  + isoRelate
  + deploidIBD 


# Bob notes
Practical: Estimating connectivity through identity measures
Simple calculation of IBS on data from two different areas - different allele frequency distributions lead to very different mean IBS
Focus on the tail. Where in space are these highly related samples? Same clusters, different clusters? Can we map the links in space?
Estimating IBD using simple ML estimator (Bob Nat Comms paper)
Estimating IBD using
hmmIBD (use OJ’s R version if easier)
isoRelate?
Plotting networks or maps
Whole genome example?
Finish with an example that ties in many of the ideas of the last 3 days
Analysis plan for a study that will sequence DR loci and neutral loci in samples from several geographic regions
Explore power to detect DR prevalence
How will you deal with polyclonal infections?
How will you look for connectivity and population structure?
Take home messages / objectives:
IBS/IBD can measure connectivity on spatial/temporal scales relevant to control
Pulling ideas from workshop together
