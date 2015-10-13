#ifndef HAPLO_CHAIN_OPTIMIZER_HPP_
#define HAPLO_CHAIN_OPTIMIZER_HPP_

class haplo_chain_optimizer {

	//Algorithm configurations
	t_count const max_iterations;

	t_position const min_overlap_length;

	bool dual_chains;

	//Model parameters
	t_haplotype haplo; //vector of size n_reads
	t_strands read_strands; //vector of size n_reads

	//Data
	alignment_data const& data;

	//Ref
	t_seq_bases const& ref;
	t_seq_bases const& alt;


	//Priors
	double ref_no_match_prior_log;
	double ref_match_prior_log;
	arma::Col<double> haplochain_log_prior;

	template<typename abort_checker>
	t_count optimize_profile(abort_checker const& ac);

	template<typename abort_checker>
	t_count optimize_pair_profile(abort_checker const& ac);

	template<typename abort_checker>
	t_count optimize_haplochains(abort_checker const& ac);

	//TODO rename : read_profile_posterior
	field<t_loglike_vector> profile_posterior(
			t_index const read_number,
			t_haplotype const& feasible_haplotypes) const;

	field<t_loglike_vector> profile_posterior(
			std::vector<t_indices> read_pairs,
			t_indices const& pairs,
			t_haplotype const& feasible_haplotypes) const;

	t_loglike_vector profile_posterior(
			t_indices const& reads,
			t_strands const& strands,
			t_haplotype const& feasible_haplotypes) const;

	//TODO rename to: find_overlapping_haplochains
	t_haplotype compute_feasible_haplotypes(
			t_position const start,
			t_position const end) const;

	double compute_haplochain_prior(t_count const unique_chains) const;

	double compute_logsum(
			t_indices const& reads,
			t_strands const& strands,
			t_position pos) const;


	double compute_posterior(
			t_haplotype const& h,
			t_strands const& strands) const;

	double compute_delta_posterior(
			t_index const read_number,
			t_strand const strand,
			t_haplochain const old_chain) const;

	double compute_delta_posterior(
			t_indices const& reads,
			t_strands const& strands,
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

	t_position read_start_postion(t_index const read_number) const;
	t_position read_end_postion(t_index const read_number) const;

	//TODO note duplicate in genotype_optimizer
	t_position haplo_chain_start(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		return data.reads_start_postions(min(reads_in_chain));
	}

	//TODO note duplicate in genotype_optimizer
	t_position haplo_chain_end(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		t_positions end_pos = data.reads_end_postions(reads_in_chain);
		return max(end_pos);
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

};

//TODO if we use logsum instead of max_logsum then we do not need ref and alt
inline haplo_chain_optimizer::haplo_chain_optimizer(
		AlgorithmConfiguration const& config,
		alignment_data const& data,
		t_seq_bases const& ref,
		t_seq_bases const& alt) :
				max_iterations(config.max_iterations),
				min_overlap_length(config.min_overlap_length),
				dual_chains(config.dual_chains),
				haplo(data.n_reads),
				read_strands(data.n_reads),
				data(data),
				ref(ref),
				alt(alt),
				haplochain_log_prior(config.haplochain_log_prior) {

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

}


template<typename abort_checker>
inline void haplo_chain_optimizer::run(const abort_checker& ac) {

	TIMER_START

	t_count change_count = 1;

	while(change_count != 0) {

		cout << "---------------------------------" << endl;

		change_count = optimize_profile(ac);

		cout << change_count << " : " << compute_posterior(haplo, read_strands) << endl;

		change_count += optimize_haplochains(ac);

		cout << change_count << " : " << compute_posterior(haplo, read_strands) << endl;


	}

}




template<typename abort_checker>
inline void haplo_chain_optimizer::run_paired(const abort_checker& ac) {

	TIMER_START

	t_count change_count = 1;

	while(change_count != 0) {

		cout << "---------------------------------" << endl;

		change_count = optimize_pair_profile(ac);

		cout << change_count << " : " << compute_posterior(haplo, read_strands) << endl;

		change_count += optimize_haplochains(ac);

		cout << change_count << " : " << compute_posterior(haplo, read_strands) << endl;

	}


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

			t_haplotype feasible_haplotypes = compute_feasible_haplotypes(read_start_postion(read_number), read_end_postion(read_number));

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

	std::vector<std::string> pair_names;
	std::vector<t_indices> read_pairs;
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

			t_position start = min(read_start_postion(read_number_1), read_start_postion(read_number_2));
			t_position end = max(read_end_postion(read_number_1), read_end_postion(read_number_2));

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


template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_haplochains(const abort_checker& ac) {

	t_count i = 0;
	for (; i < max_iterations; i++) {

		if (ac.is_aborted()) {
			return 0;
		}

		t_count changes = 0;

		t_haplotype unique_chains = unique(haplo);

		for (t_index j = 0; j < unique_chains.n_elem; ++j) {

			//Check for ctrl C
			if (ac.check_abort()) {
				break;
			}

			t_haplochain chain = unique_chains(j);

			//Find chain start and end position
			t_indices reads_in_chain = find(haplo == chain);

			t_position chain_start = read_start_postion(reads_in_chain(0));
			t_position chain_end = max(data.reads_end_postions(reads_in_chain));

			t_haplotype feasible_haplotypes = compute_feasible_haplotypes(chain_start, chain_end);

			if(feasible_haplotypes.is_empty()) {
				continue;
			}

			t_loglike_vector posterior = profile_posterior(reads_in_chain, read_strands(reads_in_chain), feasible_haplotypes);

			if (!is_finite(posterior)) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i = argmax(posterior);

			//check if any improvement
			if (1e-5 >= posterior(i)) { //TODO configable threshold
				continue;
			}

			//Update haplotype chain of read
			haplo(reads_in_chain).fill(feasible_haplotypes(i));

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

t_loglike_vector haplo_chain_optimizer::profile_posterior(
		t_indices const& reads_in_chain,
		t_strands const& strands,
		t_haplotype const& feasible_haplotypes) const {

	t_loglike_vector loglike(feasible_haplotypes.n_elem, arma::fill::zeros);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(i) = compute_delta_posterior(reads_in_chain, strands, feasible_haplotypes(i));
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

	double delta_posterior = 0;

	//Compute posterior
	new_haplo(read_number) = new_chain;
	new_strands(read_number) = strand;

	update_haplotype_chain(new_haplo, current_chain);

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

	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		delta_posterior += compute_prior_chain(new_haplo, chain);
		delta_posterior += compute_chain_loglike(new_haplo, new_strands, chain);
	}

	//TODO debug guards
//	if(fabs(delta_posterior - compute_posterior(new_haplo, new_strands) + compute_posterior(haplo, read_strands)) > 1e-5) {
//
//		cout << "chain : " << current_chain << " -> " << new_chain << endl;
//		cout << "strand : " << read_strands(read_number) << " -> " << strand << endl;
//
//		cout << delta_posterior << " vs "  << compute_posterior(new_haplo, new_strands) - compute_posterior(haplo, read_strands) << endl;
//
//		for (t_index i = 0; i < max(new_haplo); ++i) {
//
//			double a= compute_prior_chain(haplo, i);
//			double b = compute_chain_loglike(haplo, read_strands, i);
//			double c= compute_prior_chain(new_haplo, i);
//			double d = compute_chain_loglike(new_haplo, new_strands, i);
//
//			if(a != c) {
//				cout << i << " prior : " << a << " vs " << c << endl;
//			}
//
//			if(b != d) {
//				cout << i << " logsum : " << b << " vs " << d << endl;
//			}
//
//
//		}
//
//
//		throw std::runtime_error("compute_delta_posterior 1: Error");
//	}

	return delta_posterior;
}

double haplo_chain_optimizer::compute_delta_posterior(
		t_indices const& reads,
		t_strands const& strands,
		t_haplochain const new_chain) const {

	TIMER_START

	//TODO this is alot of copying
	t_haplotype new_haplo(haplo);
	t_strands new_strands(read_strands);

	t_haplochain current_chain = haplo(reads(0));

	if(current_chain == new_chain) {
		return 0;
	}

	new_haplo(reads).fill(new_chain);
	new_strands(reads) = strands;

	update_haplotype_chain(new_haplo, current_chain);

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
//	if(fabs(delta_posterior - compute_posterior(new_haplo, new_strands) + compute_posterior(haplo, read_strands)) > 1e-5) {
//
//		throw std::runtime_error("compute_delta_posterior 2: Error");
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

	if(current_chain == new_chain) {
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
	t_positions starts = data.reads_start_postions(reads_in_chain);
	t_positions ends = data.reads_end_postions(reads_in_chain);

	t_position chain_start = min(starts);
	t_position chain_end = max(ends);

	double logsum = 0;

	for (t_position pos = chain_start+1; pos <= chain_end-1; ++pos) {

		t_indices reads = reads_in_chain(find(starts <= pos-1 && pos+1 <= ends));
		logsum += compute_logsum(reads, strands, pos);

	}

	if(!std::isfinite(logsum)) {
		throw std::runtime_error("compute_chain_loglike : Error");
	}

	return logsum;
}

double haplo_chain_optimizer::compute_logsum(
		t_indices const& reads,
		t_strands const& strands,
		t_position pos) const {

	TIMER_START
	DEBUG_ENTER

	mat::fixed<6,6> logsum;

	//Prior
	for (t_epi_base a = 0; a < 6; ++a) {
		for (t_epi_base b = 0; b < 6; ++b) {

			if(a == 0 && b == 5) {
				continue;
			}

			if(a == 4 && b == 1) {
				continue;
			}

			//Ref prior
			logsum(a,b) = a % 4 + 1 == ref(pos) || a % 4 + 1 == alt(pos) ? ref_match_prior_log : ref_no_match_prior_log;
			logsum(a,b) = b % 4 + 1 == ref(pos+1) || b % 4 + 1 == alt(pos+1) ? ref_match_prior_log : ref_no_match_prior_log;
		}
	}

	//logsum


	mat::fixed<2,6> loglike_term;
	loglike_term.zeros();

	t_indices::const_iterator r = reads.begin();
	for (; r != reads.end(); ++r) {
		t_position const read_start = read_start_postion(*r);
		loglike_term += data.loglike_terms(*r, strands(*r)).rows(pos  - read_start, pos + 1 - read_start);
	}

	for (t_epi_base a = 0; a < 6; ++a) {
		logsum.row(a) += loglike_term.row(1);
	}

	for (t_epi_base b = 0; b < 6; ++b) {
		logsum.col(b) += trans(loglike_term.row(0));
	}

	logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
	logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	return logsum.max();
}



inline double haplo_chain_optimizer::compute_haplochain_prior(
		t_count const n) const {

	if(n == 0) {
		return 0;
	}

	if (n >= haplochain_log_prior.n_elem) {
		throw std::runtime_error("compute_haplochain_prior : invalid input");
	}

	return haplochain_log_prior(n-1);
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

////Find overlapping haplo chains
//inline t_haplotype haplo_chain_optimizer::compute_feasible_haplotypes(t_position const start, t_position const end) const {
//
//	t_index const read_number_start = min(data.read_numbers(start));
//	t_index const read_number_end = max(data.read_numbers(end));
//
//	t_haplotype feasible_haplotypes;
//
//	for (t_index read = read_number_start; read <= read_number_end; read++) {
//
//		//Check that read overlap haplo chain
//		if (start <= read_start_postion(read)) {
//
//			t_position overlap_length = 1 + min(end, read_end_postion(read)) - max(start, read_start_postion(read));
//
//			//TODO this is a debug check -- should not be active in release versions
//			if(overlap_length < 1) {
//				//TODO remove
//				//cout << overlap_length << endl;
//				//cout << read_end << " : " << read_end_postion(read) << " : " << read_start << " : " << read_start_postion(read) << endl;
//				throw std::runtime_error("Internal error");
//			}
//
//			if(overlap_length >= min_overlap_length) {
//				feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
//				feasible_haplotypes(feasible_haplotypes.n_elem - 1) = haplo(read);
//			}
//		}
//
//	}
//
//	return unique(feasible_haplotypes);
//}

void haplo_chain_optimizer::update_haplotype_chain(t_haplotype & haplotype, t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_indices reads_in_chain = find(haplotype == old_chain);

	t_haplochain chain = old_chain;
	t_position chain_end = read_end_postion(reads_in_chain(0));

	for (t_index i = 1; i < reads_in_chain.n_elem; ++i) {

		t_index read = reads_in_chain(i);

		if (read_start_postion(read) > chain_end - min_overlap_length + 1) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;
		}

		haplotype(read) = chain;
		chain_end = read_end_postion(read);
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
	t_position chain_end = max(read_end_postion(read_pairs[0](0)), read_end_postion(read_pairs[0](1)));

	for (t_index i = 1; i < read_pairs.size(); ++i) {

		//check if read pair is in chain
		if(haplotype(read_pairs[i](0)) != old_chain) {
			continue;
		}

		t_position start = min(read_start_postion(read_pairs[i](0)), read_start_postion(read_pairs[i](1)));

		if (start > chain_end - min_overlap_length + 1) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;
		}

		haplotype(read_pairs[i](0)) = chain;
		haplotype(read_pairs[i](1)) = chain;

		chain_end = max(read_end_postion(read_pairs[i](0)), read_end_postion(read_pairs[i](1)));
	}
}

double haplo_chain_optimizer::compute_prior_chain(t_haplotype const& h, t_haplochain chain) const {

	TIMER_START

	double reads_in_chain = sum(h == chain);

	return compute_haplochain_prior(reads_in_chain);

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


//inline t_indices haplo_chain_optimizer::get_overlapping_reads_in_chain(t_haplotype const& h,
//		t_position const pos, t_haplochain const chain) const {
//
//	TIMER_START
//
//	t_indices reads_at_pos = data.read_numbers(pos);
//	arma::uvec chain_indicator = (h(reads_at_pos) == chain);
//
//	if (sum(chain_indicator) == 0) {
//		//Return empty vector
//		return t_indices();
//	}
//
//	//Return read numbers of reads in chain
//	return reads_at_pos(find(chain_indicator));
//}

inline t_position haplo_chain_optimizer::read_start_postion(
		t_index const read_number) const {
	return data.reads_start_postions(read_number);
}

inline t_position haplo_chain_optimizer::read_end_postion(
		t_index const read_number) const {
	return data.reads_end_postions(read_number);
}

#endif /* HAPLO_CHAIN_OPTIMIZER_HPP_ */
