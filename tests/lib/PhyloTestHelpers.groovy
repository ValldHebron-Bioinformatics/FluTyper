// Shared helpers of the phylogenetics pipeline tests (tests/phylogeny.nf.test); nf-test puts tests/lib on the classpath.
class PhyloTestHelpers {

    // Only tree.nwk, phylo.log and the HTML report are published per tree; every other phylogenetics file (manifest,
    // annotations, alignment, ...) stays in the task work dir, so these files are read from the test's work dir
    // ("${launchDir}/work"; don't parse it from stdout, whose format differs between machines): the task output that
    // sits in a folder named after the tree (e.g. HA/run_manifest.json).
    // `marker` tells apart same-named files of two subtypes (e.g. H3N2 and H1N1pdm09 both have an HA tree).
    static File workFile(def workDir, String treeFolder, String name, String marker = null) {
        def root = new File(workDir.toString())
        assert root.isDirectory() : "work dir ${root} does not exist"
        def hits = []
        root.eachFileRecurse { f ->
            if (f.isFile() && f.name == name && f.parentFile.name == treeFolder && (marker == null || f.text.contains(marker))) hits << f
        }
        assert hits : "no ${treeFolder}/${name} in the work dir"
        return hits[0]
    }

    // A built tree's published folder holds exactly tree.nwk, phylo.log and the HTML report; returns phylo.log's text.
    static String checkTreeFolder(def treeDir, String treeId) {
        assert treeDir.exists() : "${treeDir} was not published"
        def names = treeDir.toFile().list().toList().sort()
        assert names == ["PhylogeneticTreeReport_${treeId}.html".toString(), "phylo.log", "tree.nwk"].sort() : "unexpected files in ${treeDir}: ${names}"
        assert treeDir.resolve("tree.nwk").text.trim().endsWith(";") : "tree.nwk is not a newick tree"
        return treeDir.resolve("phylo.log").text
    }

    // The phylogenetics steps no longer print their per-task chatter; only the summary line and actionable warnings remain
    static void checkQuietTerminal(List stdout) {
        def text = stdout.join("\n")
        ["PhyloPreflight:", "PhyloCoverageQC:", "PhyloAlignment:", "PhyloSegmentAlignment:", "PhyloTree:", "PhyloPostTree:", "PhyloReport:", "PhyloSubclades:"].each { chatter ->
            assert !text.contains(chatter) : "'${chatter}' step chatter is still printed to the terminal"
        }
    }
    static int countLines(List stdout, String needle) {
        return stdout.count { it.contains(needle) }
    }
}
