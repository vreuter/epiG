#ifndef HAPLO_CHAIN_OPTIMIZER_HPP_
#define HAPLO_CHAIN_OPTIMIZER_HPP_

class haplo_chain_optimizer {

	//Algorithm configurations
	t_count const max_iterations;

	t_position const min_overlap_length;

	//Model parameters
	t_haplotype haplo; //vector of size n_reads
	t_strands read_strands; //vector of size n_reads

	//Data
	alignment_data const& data;

	std::vector<t_indices> read_pairs;

	//Ref
	t_seq_bases const& ref;
	t_seq_bases const& alt;


	//Priors
	double ref_no_match_prior_log;
	double ref_match_prior_log;

	arma::Col<double> haplochain_log_prior;

	field<vec> ref_priores;

	template<typename abort_checker>
	t_count optimize_profile(abort_checker const& ac);

	template<typename abort_checker>
	t_count optimize_pair_profile(abort_checker const& ac);

	//TODO rename : read_profile_posterior
	field<t_loglike_vector> profile_posterior(
			t_index const read_number,
			t_haplotype const& feasible_haplotypes) const;

	field<t_loglike_vector> profile_posterior(
			std::vector<t_indices> read_pairs,
			t_indices const& pairs,
			t_haplotype const& feasible_haplotypes) const;

	//TODO rename to: find_overlapping_haplochains
	t_haplotype compute_feasible_haplotypes(
			t_position const start,
			t_position const end) const;

	field<vec> compute_ref_priors() const;

	double compute_logsum(
			mat const& loglike_term,
			t_position pos) const;


	double compute_posterior(
			t_haplotype const& h,
			t_strands const& strands) const;

	double compute_delta_posterior(
			t_index const read_number,
			t_strand const strand,
			t_haplochain const old_chain) const;

	double compute_delta_posterior(
			std::vector<t_indices> read_pairs,
			t_indices const& reads,
			t_strand const strand,
			t_haplochain const old_chain) const;

	double compute_chain_loglike(
			t_haplotype const& h,
			t_strands const& strands,
			t_haplochain const chain) const;

	double compute_chain_loglike(
			t_haplotype const& h,
			t_strands const& strands,
			t_haplochain const chain,
			t_indices const& reads_overlapping_region,
			t_position region_start,
			t_position region_end) const;

	double compute_prior_chain(
			t_haplotype const& h,
			t_haplochain chain) const;

	//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
	void update_haplotype_chain(
			t_haplotype & haplotype,
			t_haplochain const old_chain) const;

	void update_haplotype_chain(
			std::vector<t_indices> read_pairs,
			t_haplotype & haplotype,
			t_haplochain const old_chain) const;

	t_genotype compute_chain_genotype(t_haplochain const chain) const;

	void compute_genotype(
			mat const& loglike_term,
			t_position pos,
			t_epi_base & g1,
			t_epi_base & g2) const;


	t_position read_start_position(t_index const read_number) const;
	t_position read_end_position(t_index const read_number) const;

	//TODO note duplicate in genotype_optimizer
	t_position haplo_chain_start(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		return data.reads_start_positions(min(reads_in_chain));
	}

	//TODO note duplicate in genotype_optimizer
	t_position haplo_chain_end(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		t_positions end_pos = data.reads_end_positions(reads_in_chain);
		return max(end_pos);
	}

	void chain_clean() {

		t_haplotype unique_chains = sort(unique(haplo));

		t_haplotype new_haplo(haplo);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			t_indices reads_in_chain(find(haplo == unique_chains(i)));
			new_haplo(reads_in_chain).fill(i);
		}

		haplo = new_haplo;
	}

public:

	haplo_chain_optimizer(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt);

	template<typename abort_checker>
	void run(abort_checker const& ac);

	template<typename abort_checker>
	void run_paired(abort_checker const& ac);

	t_haplotype get_haplotype_chains() const {
		return haplo;
	}

	t_strands get_read_strands() const {
		return read_strands;
	}

	field<t_genotype> get_chain_genotypes() const;

	t_positions haplo_chain_start() const {

		t_haplotype unique_chains = sort(unique(haplo));

		t_positions pos(unique_chains.n_elem);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			pos(i) = haplo_chain_start(unique_chains(i));
		}

		return pos;
	}

	t_positions haplo_chain_end() const {

		t_haplotype unique_chains = sort(unique(haplo));

		t_positions pos(unique_chains.n_elem);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			pos(i) = haplo_chain_end(unique_chains(i));
		}

		return pos;
	}


};

//TODO if we use logsum instead of max_logsum then we do not need ref and alt
inline haplo_chain_optimizer::haplo_chain_optimizer(
		AlgorithmConfiguration const& config,
		alignment_data const& data,
		t_seq_bases const& ref,
		t_seq_bases const& alt) :
				max_iterations(config.max_iterations),
				min_overlap_length(config.min_overlap_length),
				haplo(data.n_reads),
				read_strands(data.n_reads),
				data(data),
				ref(ref),
				alt(alt),
				haplochain_log_prior(config.haplochain_log_prior),
				ref_priores() {

	///////////////// Initialize:

	read_strands.zeros();

	///// Init haplo chains


	if(config.use_paired_reads) {

		// assign each read pair to a chain
		haplo.zeros();


		std::vector<std::string> pair_names;
		t_index i = 1;
		for(t_index read_number = 0; read_number < data.n_reads; ++read_number) {

			std::vector<std::string>::const_iterator it = find(data.read_names.begin()+read_number+1, data.read_names.end(), data.read_names[read_number]);

			if(it != data.read_names.end()) {
				haplo(read_number) = i;
				haplo(it - data.read_names.begin()) = i;
				++i;
			}

			else {

				if(haplo(read_number) == 0) {
					haplo(read_number) = i;
					++i;
				}
			}


		}

		haplo = haplo -1;

	} else {
		//assign each read to a chain
		for(t_index i = 0; i < haplo.n_elem; ++i) {
			haplo(i) = i;
		}
	}

	//Init ref prior
	ref_no_match_prior_log = log(1-config.ref_prior);
	ref_match_prior_log = log(config.ref_prior);

	ref_priores = compute_ref_priors();

	//Find read pairs
	// and init read_pairs

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

}


template<typename abort_checker>
inline void haplo_chain_optimizer::run(const abort_checker& ac) {

	TIMER_START

	cout << "---------------------------------" << endl;

	t_count change_count = optimize_profile(ac);

	cout << change_count << " : " << compute_posterior(haplo, read_strands) << endl;

	chain_clean();

}




template<typename abort_checker>
inline void haplo_chain_optimizer::run_paired(const abort_checker& ac) {

	TIMER_START

	cout << "---------------------------------" << endl;

	t_count change_count = optimize_pair_profile(ac);

	cout << change_count << " : " << compute_posterior(haplo, read_strands) << endl;

	chain_clean();

}



template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_profile(const abort_checker& ac) {

	t_count i = 0;
	for (; i < max_iterations; i++) {

		if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		for (t_index read_number = 0; read_number < data.n_reads;
				++read_number) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			t_haplotype feasible_haplotypes = compute_feasible_haplotypes(read_start_position(read_number), read_end_position(read_number));

			//Add free haplo chain , i.e. a haplo chain with no reads
			feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
			feasible_haplotypes(feasible_haplotypes.n_elem - 1) = max(haplo) + 1;

			//field 0 - fwd , 1 - rev
			field<t_loglike_vector> loglike = profile_posterior(read_number, feasible_haplotypes);

			if (!is_finite(loglike(0)) || !is_finite(loglike(1))) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i_fwd = argmax(loglike(0));
			t_index i_rev = argmax(loglike(1));

			//check if any improvement
			if (1e-5 >= loglike(0)(i_fwd) && 1e-5 >= loglike(1)(i_rev)) { //TODO configable threshold
				continue;
			}

			//Get current haplotype chain of read
			t_haplochain old_chain = haplo(read_number);

			//Update haplotype chain of read
			if(loglike(0)(i_fwd) > loglike(1)(i_rev)) {
				haplo(read_number) = feasible_haplotypes(i_fwd);
				read_strands(read_number) = strand_fwd;
			}

			else {
				haplo(read_number) = feasible_haplotypes(i_rev);
				read_strands(read_number) = strand_rev;
			}

			//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
			update_haplotype_chain(haplo, old_chain);

			//add change
			changes++;

		}

		//TODO check posterior is decreasing
		//TODO remove output
		//cout << changes << " : " << compute_posterior(haplo) << endl;

		if (changes == 0) {
			break;
		}

	}

	if (i == max_iterations) {
		report_error("Max iteration limit reached");
	}

	return i;

}

template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_pair_profile(
		const abort_checker& ac) {

	t_count i = 0;
	for (; i < max_iterations; i++) {

		if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		for (t_index pair_number = 0; pair_number < read_pairs.size(); ++pair_number) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			t_haplotype const& reads_pair = read_pairs[pair_number];

			t_index read_number_1 = reads_pair(0);
			t_index read_number_2 = reads_pair(1);

			//TODO debug guards
			if(haplo(read_number_1) != haplo(read_number_2)) {
				throw std::runtime_error("optimize_pair_profile : internal error : haplotype chain structure corrupt");
			}

			t_position start = min(read_start_position(read_number_1), read_start_position(read_number_2));
			t_position end = max(read_end_position(read_number_1), read_end_position(read_number_2));

			t_haplotype feasible_haplotypes = compute_feasible_haplotypes(start, end);

			//Add free haplo chain , i.e. a haplo chain with no reads
			feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
			feasible_haplotypes(feasible_haplotypes.n_elem - 1) = max(haplo) + 1;

			//field 0 - fwd , 1 - rev
			field<t_loglike_vector> loglike = profile_posterior(read_pairs, reads_pair, feasible_haplotypes);


			if (!is_finite(loglike(0)) || !is_finite(loglike(1))) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i_fwd = argmax(loglike(0));
			t_index i_rev = argmax(loglike(1));

			//check if any improvement
			if (1e-5 >= loglike(0)(i_fwd) && 1e-5 >= loglike(1)(i_rev)) { //TODO configable threshold
				continue;
			}


			//Get current haplotype chain of read pair
			t_haplochain old_chain = haplo(read_number_1);

			//Update haplotype chain of read pair
			if(loglike(0)(i_fwd) > loglike(1)(i_rev)) {
				haplo(read_number_1) = feasible_haplotypes(i_fwd);
				haplo(read_number_2) = feasible_haplotypes(i_fwd);
				read_strands(read_number_1) = strand_fwd;
				read_strands(read_number_2) = strand_fwd;
			}

			else {
				haplo(read_number_1) = feasible_haplotypes(i_rev);
				haplo(read_number_2) = feasible_haplotypes(i_rev);
				read_strands(read_number_1) = strand_rev;
				read_strands(read_number_2) = strand_rev;
			}

			//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
			update_haplotype_chain(read_pairs, haplo, old_chain);

			//add change
			changes++;

		}

		if (changes == 0) {
			break;
		}
	}

	if (i == max_iterations) {
		report_error("Max iteration limit reached");
	}

	return i;

}

field<t_loglike_vector> haplo_chain_optimizer::profile_posterior(
		t_index const read_number,
		t_haplotype const& feasible_haplotypes) const {

	DEBUG_ENTER

	field<t_loglike_vector> loglike(2);
	loglike(0).zeros(feasible_haplotypes.n_elem);
	loglike(1).zeros(feasible_haplotypes.n_elem);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(0)(i) = compute_delta_posterior(read_number, strand_fwd, feasible_haplotypes(i));
		loglike(1)(i) = compute_delta_posterior(read_number, strand_rev, feasible_haplotypes(i));
	}

	return loglike;
}


field<t_loglike_vector> haplo_chain_optimizer::profile_posterior(
		std::vector<t_indices> read_pairs,
		t_indices const& pair,
		t_haplotype const& feasible_haplotypes) const {

	field<t_loglike_vector> loglike(2);
	loglike(0).zeros(feasible_haplotypes.n_elem);
	loglike(1).zeros(feasible_haplotypes.n_elem);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(0)(i)  = compute_delta_posterior(read_pairs, pair, strand_fwd, feasible_haplotypes(i));
		loglike(1)(i)  = compute_delta_posterior(read_pairs, pair, strand_rev, feasible_haplotypes(i));
	}

	return loglike;
}

double haplo_chain_optimizer::compute_delta_posterior(
		t_index const read_number,
		t_strand strand,
		t_haplochain const new_chain) const {

	TIMER_START
	DEBUG_ENTER

	t_haplotype new_haplo(haplo); //TODO this is alot of copying
	t_strands new_strands(read_strands); //TODO this is alot of copying

	t_haplochain current_chain = haplo(read_number);

	if(current_chain == new_chain && strand == read_strands(read_number)) {
		return 0;
	}

	t_position region_start = read_start_position(read_number)-1;
	t_position region_end = read_end_position(read_number)+1;

	t_indices reads_in_region = find(data.reads_start_positions <= region_end && data.reads_end_positions >= region_start);

	double delta_posterior = 0;

	//Compute posterior
	new_haplo(read_number) = new_chain;
	new_strands(read_number) = strand;

	update_haplotype_chain(new_haplo, current_chain);

	//current posterior
	delta_posterior -= compute_prior_chain(haplo, current_chain);
	delta_posterior -= compute_chain_loglike(haplo, read_strands, current_chain, reads_in_region, region_start, region_end);

	//new posterior
	delta_posterior += compute_prior_chain(new_haplo, current_chain);
	delta_posterior += compute_chain_loglike(new_haplo, new_strands, current_chain, reads_in_region, region_start, region_end);

	if(new_chain != current_chain) {
		//current posterior
		delta_posterior -= compute_prior_chain(haplo, new_chain);
		delta_posterior -= compute_chain_loglike(haplo, read_strands, new_chain, reads_in_region, region_start, region_end);

		//new posterior
		delta_posterior += compute_prior_chain(new_haplo, new_chain);
		delta_posterior += compute_chain_loglike(new_haplo, new_strands, new_chain, reads_in_region, region_start, region_end);

	}

	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		delta_posterior += compute_prior_chain(new_haplo, chain);
		delta_posterior += compute_chain_loglike(new_haplo, new_strands, chain, reads_in_region, region_start, region_end);

	}

	//TODO debug guards
//	if(fabs(delta_posterior - compute_posterior(new_haplo, new_strands) + compute_posterior(haplo, read_strands)) > 1e-5) {
//
//		cout << "chain : " << current_chain << " -> " << new_chain << endl;
//		cout << "strand : " << read_strands(read_number) << " -> " << strand << endl;
//
//		cout << delta_posterior << " vs "  << compute_posterior(new_haplo, new_strands) - compute_posterior(haplo, read_strands) << endl;
//
//		for (t_index i = 0; i <= max(new_haplo); ++i) {
//
//			double a = compute_prior_chain(haplo, i);
//			double b = compute_chain_loglike(haplo, read_strands, i);
//			double c = compute_prior_chain(new_haplo, i);
//			double d = compute_chain_loglike(new_haplo, new_strands, i);
//			double e = compute_chain_loglike(haplo, read_strands, i, reads_in_region, region_start, region_end);
//			double f = compute_chain_loglike(new_haplo, new_strands, i, reads_in_region, region_start, region_end);
//
//			if(a != c) {
//				cout << i << " prior : " << a << " vs " << c << endl;
//			}
//
//			if(b != d || e != f ) {
//				cout << " ---------------------------- " << endl;
//				cout << i << " logsum : " << b << " vs " << d << endl;
//				cout << i << " logsum : " << e << " vs " << f << endl;
//
//			}
//
//
//		}
//
//		throw std::runtime_error("compute_delta_posterior 1: Error");
//	}

	return delta_posterior;
}

double haplo_chain_optimizer::compute_delta_posterior(
		std::vector<t_indices> read_pairs,
		t_indices const& pair,
		t_strand const strand,
		t_haplochain const new_chain) const {

	TIMER_START

	t_haplotype new_haplo(haplo);
	t_strands new_strands(read_strands); //TODO this is alot of copying

	t_haplochain current_chain = haplo(pair(0));

	if(current_chain == new_chain && strand == read_strands(pair(0))) {
		return 0;
	}

	//new posterior
	new_haplo(pair).fill(new_chain);
	new_strands(pair).fill(strand);
	update_haplotype_chain(read_pairs, new_haplo, current_chain);

	double delta_posterior = 0;

	//current posterior
	delta_posterior -= compute_prior_chain(haplo, current_chain);
	delta_posterior -= compute_chain_loglike(haplo, read_strands, current_chain);

	//new posterior
	delta_posterior += compute_prior_chain(new_haplo, current_chain);
	delta_posterior += compute_chain_loglike(new_haplo, new_strands, current_chain);

	if(new_chain != current_chain) {
		//current posterior
		delta_posterior -= compute_prior_chain(haplo, new_chain);
		delta_posterior -= compute_chain_loglike(haplo, read_strands, new_chain);

		//new posterior
		delta_posterior += compute_prior_chain(new_haplo, new_chain);
		delta_posterior += compute_chain_loglike(new_haplo, new_strands, new_chain);
	}

	//add likelihood for new chains
	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		delta_posterior += compute_prior_chain(new_haplo, chain);
		delta_posterior += compute_chain_loglike(new_haplo, new_strands, chain);
	}

	//TODO debug guards
//		if(fabs(delta_posterior - compute_posterior(new_haplo, new_strands) + compute_posterior(haplo, read_strands)) > 1e-5) {
//
//			throw std::runtime_error("compute_delta_posterior 3: Error");
//		}


	return delta_posterior;
}

double haplo_chain_optimizer::compute_chain_loglike(
		t_haplotype const& h,
		t_strands const& strands,
		t_haplochain const chain) const {

	TIMER_START
	DEBUG_ENTER

	t_indices reads_in_chain = find(h == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	//find start, end pos of chain
	t_positions starts = data.reads_start_positions(reads_in_chain);
	t_positions ends = data.reads_end_positions(reads_in_chain);

	t_position chain_start = min(starts);
	t_position chain_end = max(ends);

	double logsum = 0;

	mat loglike_term(chain_end - chain_start + 1, 6, fill::zeros);

	t_indices::const_iterator r = reads_in_chain.begin();
	for (; r != reads_in_chain.end(); ++r) {
		t_position const read_start = read_start_position(*r);
		t_position const read_end = read_end_position(*r);

		loglike_term.rows(read_start - chain_start, read_end - chain_start)
				+= data.loglike_terms(*r, strands(*r));
	}

	for (t_position pos = chain_start; pos <= chain_end-1; ++pos) {
		logsum += compute_logsum(loglike_term.rows(pos - chain_start, pos - chain_start + 1), pos);
	}

	if(!std::isfinite(logsum)) {
		throw std::runtime_error("compute_chain_loglike : Error");
	}

	return logsum;
}

double haplo_chain_optimizer::compute_chain_loglike(
		t_haplotype const& h,
		t_strands const& strands,
		t_haplochain const chain,
		t_indices const& reads_overlapping_region,
		t_position region_start,
		t_position region_end) const {

	TIMER_START
	DEBUG_ENTER

	if (reads_overlapping_region.n_elem == 0) {
		return 0;
	}

	t_indices reads_in_chain = find(h(reads_overlapping_region) == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	reads_in_chain = reads_overlapping_region(reads_in_chain);

	//find start, end pos of chain
	t_positions starts = data.reads_start_positions(reads_in_chain);
	t_positions ends = data.reads_end_positions(reads_in_chain);

	t_position chain_start = min(starts);
	t_position chain_end = max(ends);

	double logsum = 0;

	mat loglike_term(chain_end - chain_start + 1, 6, fill::zeros);

	t_indices::const_iterator r = reads_in_chain.begin();
	for (; r != reads_in_chain.end(); ++r) {
		t_position const read_start = read_start_position(*r);
		t_position const read_end = read_end_position(*r);

		loglike_term.rows(read_start - chain_start, read_end - chain_start)
				+= data.loglike_terms(*r, strands(*r));
	}

	if(region_end >= chain_end) {
		region_end = chain_end - 1;
	}

	region_start -= 1;
	if(region_start < chain_start) {
		region_start = chain_start;
	}

	for (t_position pos = region_start; pos <= region_end; ++pos) {
		logsum += compute_logsum(loglike_term.rows(pos - chain_start, pos - chain_start + 1), pos);
	}

	if(!std::isfinite(logsum)) {
		throw std::runtime_error("compute_chain_loglike : Error");
	}

	return logsum;
}


field<vec> haplo_chain_optimizer::compute_ref_priors() const {

	field<vec> priores(5,5);

	for (t_base r = 0; r < 5; ++r) {
		for (t_base a = 0; a < 5; ++a) {
			priores(r,a).set_size(6);
			for (t_epi_base g = 0; g < 6; ++g) {
				priores(r, a)(g) = (g % 4 + 1 == r || g % 4 + 1 == a ? ref_match_prior_log : ref_no_match_prior_log);
			}
		}
	}


	return priores;
}

double haplo_chain_optimizer::compute_logsum(
		mat const& loglike_term,
		t_position pos) const {

	DEBUG_ENTER
	//TIMER_START

	//Prior scaled with number of overlapping reads
	//Logsum

	mat::fixed<6,6> logsum;

	logsum.each_row() = trans(ref_priores(ref(pos+1), alt(pos+1))) + loglike_term.row(1);
	logsum.each_col() += ref_priores(ref(pos), alt(pos)) + trans(loglike_term.row(0));

	//rule out half methylated genotypes
	logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
	logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	return logsum.max();
}

t_genotype haplo_chain_optimizer::compute_chain_genotype(t_haplochain const chain) const {

	t_indices reads_in_chain = find(haplo == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	//find start, end pos of chain
	t_positions starts = data.reads_start_positions(reads_in_chain);
	t_positions ends = data.reads_end_positions(reads_in_chain);

	t_position chain_start = min(starts);
	t_position chain_end = max(ends);

	t_genotype genotype(chain_end-chain_start+1, fill::zeros);

	mat loglike_term(chain_end - chain_start + 1, 6, fill::zeros);

	t_indices::const_iterator r = reads_in_chain.begin();
	for (; r != reads_in_chain.end(); ++r) {
		t_position const read_start = read_start_position(*r);
		t_position const read_end = read_end_position(*r);

		loglike_term.rows(read_start - chain_start, read_end - chain_start) += data.loglike_terms(*r, read_strands(*r));
	}

	t_epi_base g;
	for (t_position pos = chain_start; pos <= chain_end-1; ++pos) {
		compute_genotype(
				loglike_term.rows(pos - chain_start, pos - chain_start + 1),
				pos,
				genotype(pos-chain_start),
				g);
	}

	genotype(chain_end - chain_start) = g;

	return genotype;

}

void haplo_chain_optimizer::compute_genotype(
		mat const& loglike_term,
		t_position pos,
		t_epi_base & g1,
		t_epi_base & g2) const {

	DEBUG_ENTER
	//TIMER_START

	//Prior scaled with number of overlapping reads
	//Logsum

	mat::fixed<6,6> logsum;

	logsum.each_row() = trans(ref_priores(ref(pos+1), alt(pos+1))) + loglike_term.row(1);
	logsum.each_col() += ref_priores(ref(pos), alt(pos)) + trans(loglike_term.row(0));

	//rule out half methylated genotypes
	logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
	logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	t_index a;
	t_index b;
	logsum.max(a, b);

	g1 = static_cast<t_epi_base>(a);
	g2 = static_cast<t_epi_base>(b);

}


field<t_genotype> haplo_chain_optimizer::get_chain_genotypes() const {

	t_haplotype unique_chains = sort(unique(haplo));

	field<t_genotype> genotypes(unique_chains.n_elem);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			genotypes(i) = compute_chain_genotype(unique_chains(i));
		}

	return(genotypes);
}

t_haplotype haplo_chain_optimizer::compute_feasible_haplotypes(t_position const start, t_position const end) const {

	TIMER_START

	t_haplotype feasible_haplotypes;

	//Compute for each chain the number bases overlapping
	t_haplotype unique_chains = unique(haplo);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {
		t_position chain_start = haplo_chain_start(unique_chains(i));
		t_position chain_end = haplo_chain_end(unique_chains(i));

		t_position overlap_length = pos(end-start + 1 - pos(chain_start - start) - pos(end-chain_end));

		if(overlap_length >= min_overlap_length) {
			feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
			feasible_haplotypes(feasible_haplotypes.n_elem - 1) = unique_chains(i);
		}

	}

	return feasible_haplotypes;

}

void haplo_chain_optimizer::update_haplotype_chain(t_haplotype & haplotype, t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_indices reads_in_chain = find(haplotype == old_chain);

	t_haplochain chain = old_chain;
	t_position chain_end = read_end_position(reads_in_chain(0));

	for (t_index i = 1; i < reads_in_chain.n_elem; ++i) {

		t_index read = reads_in_chain(i);

		t_position overlap_length = chain_end - read_start_position(read) + 1;
		if(overlap_length < min_overlap_length) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;

			haplotype(read) = chain;
			chain_end = read_end_position(read);
		}

		else {
			haplotype(read) = chain;
			chain_end = max(chain_end, read_end_position(read));
		}
	}
}

void haplo_chain_optimizer::update_haplotype_chain(
		std::vector<t_indices> read_pairs,
		t_haplotype & haplotype,
		t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_haplochain chain = old_chain;
	t_position chain_end = max(read_end_position(read_pairs[0](0)), read_end_position(read_pairs[0](1)));

	for (t_index i = 1; i < read_pairs.size(); ++i) {

		//check if read pair is in chain
		if(haplotype(read_pairs[i](0)) != old_chain) {
			continue;
		}

		t_position start = min(read_start_position(read_pairs[i](0)), read_start_position(read_pairs[i](1)));

		if (start > chain_end - min_overlap_length + 1) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;

			haplotype(read_pairs[i](0)) = chain;
			haplotype(read_pairs[i](1)) = chain;

			chain_end = max(read_end_position(read_pairs[i](0)), read_end_position(read_pairs[i](1)));

		}

		else {
			haplotype(read_pairs[i](0)) = chain;
			haplotype(read_pairs[i](1)) = chain;

			t_position read_pair_end = max(read_end_position(read_pairs[i](0)), read_end_position(read_pairs[i](1)));
			chain_end = max(chain_end, read_pair_end);
		}
	}
}

double haplo_chain_optimizer::compute_prior_chain(t_haplotype const& h, t_haplochain chain) const {

	TIMER_START

	t_count reads_in_chain = sum(h == chain);

	if(reads_in_chain == 0) {
		return 0;
	}

	return haplochain_log_prior(reads_in_chain-1);

}

double haplo_chain_optimizer::compute_posterior(
		t_haplotype const& h,
		t_strands const& strands) const {

	TIMER_START

	//Compute prior

	double posterior = 0;

	// Compute likelihood
	t_haplotype unique_chains = unique(h);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {
		posterior += compute_prior_chain(h, unique_chains(i));
		posterior += compute_chain_loglike(h, strands, unique_chains(i));

	}

	return posterior;
}

inline t_position haplo_chain_optimizer::read_start_position(
		t_index const read_number) const {
	return data.reads_start_positions(read_number);
}

inline t_position haplo_chain_optimizer::read_end_position(
		t_index const read_number) const {
	return data.reads_end_positions(read_number);
}

#endif /* HAPLO_CHAIN_OPTIMIZER_HPP_ */
