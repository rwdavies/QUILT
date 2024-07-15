* v2.0.0
	* Major version update
	* Includes msPBWT method of haplotype selection
	* Includes common and all SNP processing
	* Includes NIPT method for imputing mother and fetus from cfDNA NIPT sequence
* v1.0.5
	* Be able to work with cram files
* v1.0.4
	* Mostly small bugs and some future capabilities
* v1.0.3
	* Simplify how QUILT HLA reference packages are built
* v1.0.2
	* Fix minor but that prevented HLA reference panel from building on new reference data
* v1.0.1
	* Fix minor but that prevented HLA reference panel from building
* v1.0.0
	* Bump major version number as first version post paper release
	* Fix bugs related to rare haplotypes in full panel HMM
* v0.1.9
	* Change default output of GT entry to phased value to facilitate ligating phase reslts together
	* Add code to try and re-run samples with different parameter values when underflow happens
* v0.1.8
	* RAM decrease when running many samples
* v0.1.7
	* Speedups and RAM decrease when using reference panels with fewer haplotypes
* v0.1.6
	* Discovered bug where only 1 read caused an error. Fix (obviate) this by requiring 2 of more reads through parameter minimum_number_of_sample_reads with default currently 2, where samples with fewer than 2 reads have all missing output
* v0.1.5
	* Export downsampleToCov and set to 30 to reduce underflow likelihood for high coverage regions
	* Add minGLValue as bound on haplotype genotype likelihoods also to reduce underflow likelihood
* v0.1.4
	* Option to impute without specifying genetic map, just using expected genetic and physical distance
* v0.1.3
	* Bugfixes
* v0.1.2
	* More optimization, avoid normalization when possible, avoid unnecessary adds and mults
* v0.1.1
	* Better acceptance tests
	* More validation
	* Use of genfile and/or phasefile
* v0.1.0
	* First versionned beta release
* v0.0.0
	* Generic development and pre-release version
