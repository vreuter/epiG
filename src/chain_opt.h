#ifndef HAPLO_CHAIN_OPTIMIZER_HPP_
#define HAPLO_CHAIN_OPTIMIZER_HPP_

class haplo_chain_optimizer {

	//Algorithm configurations
	AlgorithmConfiguration const& config;
	bool use_paired_reads;

	//Data
	alignment_data const& data;
	t_seq_bases const& ref;
	t_seq_bases const& alt;

	//Holds haplotype model
	haplotype h;

	//max iterations
	t_count const max_iterations;

	//Report warnings to R
	omp_rwarn & warnings;

	template<typename abort_checker>
	t_count optimize_profile(abort_checker const& ac);

	field<t_loglike_vector> profile_posterior(
			t_index const id,
			t_haplochains const& feasible_haplotypes) const;

	std::vector<t_indices> find_read_pairs(alignment_data const& data) const;

public:

	haplo_chain_optimizer(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			omp_rwarn & warnings,
			bool use_paired_reads_dummy);

	haplo_chain_optimizer(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			omp_rwarn & warnings);

	template<typename abort_checker>
	void run(abort_checker const& ac);

	t_haplochains get_haplochains() const {
		return h.haplotypes();
	}

	t_strands get_read_strands() const  {
		return h.strands();
	}

	field<t_genotype> get_chain_genotypes() const {
		return h.genotypes();
	}

	field<mat> get_chain_loglikes() const {
			return h.loglikes();
	}


	t_positions haplo_chain_start() const {
		return h.starts();
	}

	t_positions haplo_chain_end() const {
		return h.ends();
	}


	static haplo_chain_optimizer create(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			omp_rwarn & warnings,
			bool use_paired_reads) {

		DEBUG_ENTER

		if(use_paired_reads) {
			return haplo_chain_optimizer(config, data, ref, alt, warnings, use_paired_reads);
		}

		return haplo_chain_optimizer(config, data, ref, alt, warnings);
	}
};


//TODO if we use logsum instead of max_logsum then we do not need ref and alt
inline haplo_chain_optimizer::haplo_chain_optimizer(
		AlgorithmConfiguration const& config,
		alignment_data const& data,
		t_seq_bases const& ref,
		t_seq_bases const& alt,
		omp_rwarn & warnings,
		bool use_paired_reads_dummy) :
				config(config),
				use_paired_reads(true),
				data(data),
				ref(ref),
				alt(alt),
				h(config, data, ref, alt, config.min_overlap_length, find_read_pairs(data)),
				max_iterations(config.max_iterations),
				warnings(warnings)
		{
}

inline haplo_chain_optimizer::haplo_chain_optimizer(
		AlgorithmConfiguration const& config,
		alignment_data const& data,
		t_seq_bases const& ref,
		t_seq_bases const& alt,
		omp_rwarn & warnings) :
				config(config),
				use_paired_reads(false),
				data(data),
				ref(ref),
				alt(alt),
				h(config, data, ref, alt, config.min_overlap_length),
				max_iterations(config.max_iterations),
				warnings(warnings)
		{
}

template<typename abort_checker>
inline void haplo_chain_optimizer::run(const abort_checker& ac) {

	DEBUG_ENTER
	TIMER_START

	std::vector<t_indices> read_pairs;

	//TODO fix double work : some of this work was already done in the init

	if(use_paired_reads) {
		read_pairs = find_read_pairs(data);
	}

	for(int stage = 1; stage < config.max_stages; stage++) {


		if(use_paired_reads) {
			h.set_blocks(read_pairs);
		}

		else {
			h.set_blocks();
		}

		optimize_profile(ac);

		h.chain_clean();
		h.set_blocks(h);
		//h.set_min_overlap(); TODO different min_overlap for stages

		t_count change_count = optimize_profile(ac);
		h.chain_clean();

		//cout << change_count << " : " << h.posterior() << endl;

		if(change_count <= 0) {
			break;
		}

	}

}

template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_profile(
		const abort_checker& ac) {

	DEBUG_ENTER

	t_count i = 0;
	for (; i < max_iterations; i++) {

		if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		for (t_index id = 0; id < h.n_elements; ++id) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			t_haplochains feasible_haplotypes = h.compute_feasible_haplotypes(id);

			//field 0 - fwd , 1 - rev
			field<t_loglike_vector> loglike = profile_posterior(id, feasible_haplotypes);

			if (!is_finite(loglike(0)) || !is_finite(loglike(1))) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i_fwd = argmax(loglike(0));
			t_index i_rev = argmax(loglike(1));

			//check if any improvement
			if (1e-7 >= loglike(0)(i_fwd) && 1e-7 >= loglike(1)(i_rev)) { //TODO configable threshold
				continue;
			}

			//Update haplotype chain
			if(loglike(0)(i_fwd) > loglike(1)(i_rev)) {

//				if(1e-5 >= loglike(0)(i_fwd)) {
//					//No improvement
//
//				}

				h.move(id, strand_fwd, feasible_haplotypes(i_fwd));
			}

			else {
				h.move(id, strand_rev, feasible_haplotypes(i_rev));
			}

			//add change
			changes++;

		}

		if (changes == 0) {
			break;
		}
	}

	if (i == max_iterations) {
		warnings.add("Max iteration limit reached");
	}

	return i;

}

field<t_loglike_vector> haplo_chain_optimizer::profile_posterior(
		t_index const id,
		t_haplochains const& feasible_haplotypes) const {

	DEBUG_ENTER

	field<t_loglike_vector> loglike(2);
	loglike(0).zeros(feasible_haplotypes.n_elem);
	loglike(1).zeros(feasible_haplotypes.n_elem);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(0)(i) = h.delta_posterior(id, strand_fwd, feasible_haplotypes(i));
		loglike(1)(i) = h.delta_posterior(id, strand_rev, feasible_haplotypes(i));
	}

	return loglike;
}

std::vector<t_indices> haplo_chain_optimizer::find_read_pairs(alignment_data const& data) const {

	DEBUG_ENTER

	std::vector<t_indices> read_pairs;
	std::vector<std::string> pair_names;

	for(t_index read_number = 0; read_number < data.n_reads; ++read_number) {

		std::vector<std::string>::iterator it = find(pair_names.begin(), pair_names.end(), data.read_names[read_number]);

		if(it != pair_names.end()) {
			t_indices & rp = read_pairs[it - pair_names.begin()];
			rp(1) = read_number;
		}

		else {
			pair_names.push_back(data.read_names[read_number]);
			t_indices read_pair(2);
			read_pair.fill(read_number);
			read_pairs.push_back(read_pair);
		}
	}

	return read_pairs;
}
#endif /* HAPLO_CHAIN_OPTIMIZER_HPP_ */
