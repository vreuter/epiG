#ifndef HAPLO_CHAIN_OPTIMIZER_HPP_
#define HAPLO_CHAIN_OPTIMIZER_HPP_

class haplo_chain_optimizer {

	//Algorithm configurations
	t_count const max_iterations;

	t_position const min_overlap_length;

	//Model parameters
	t_haplotype haplo; //vector of size n_reads
	t_strands strands;

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
	t_loglike_vector profile_posterior(t_index const read_number,
			t_haplotype const& feasible_haplotypes) const;

	t_loglike_vector profile_posterior(t_indices const& reads, t_haplotype const& feasible_haplotypes) const;

	t_loglike_vector profile_posterior(std::vector<t_indices> read_pairs, t_indices const& pair, t_haplotype const& feasible_haplotypes) const;

	//TODO rename to: find_overlapping_haplochains
	t_haplotype compute_feasible_haplotypes(t_position const start, t_position const end) const;

	double compute_haplochain_prior(t_count const unique_chains) const;

	t_loglike_vector compute_logsum(t_indices const& reads, t_position const pos) const;

	double compute_posterior(t_haplotype const& h) const;

	double compute_delta_posterior(t_index const read_number, t_haplochain const old_chain) const;
	double compute_delta_posterior(t_indices const& reads, t_haplochain const old_chain) const;
	double compute_delta_posterior(std::vector<t_indices> read_pairs, t_indices const& reads, t_haplochain const old_chain) const;

	double compute_chain_loglike(t_haplotype const& h, t_haplochain const chain) const;
	//double compute_prior(t_haplotype const& h, t_position const start, t_position const end) const;
	double compute_prior_chain(t_haplotype const& h, t_haplochain chain) const;

	//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
	void update_haplotype_chain(t_haplotype & haplotype, t_haplochain const old_chain) const;
	void update_haplotype_chain(std::vector<t_indices> read_pairs, t_haplotype & haplotype, t_haplochain const old_chain) const;

	t_strand chain_strand(t_haplochain const chain) const;

	t_indices get_overlapping_reads_in_chain(t_haplotype const& h, t_position const pos,
			t_haplochain const chain) const;

	t_position read_start_postion(t_index const read_number) const;
	t_position read_end_postion(t_index const read_number) const;

	//TODO note duplicate in genotype_optimizer
	t_position haplo_chain_start(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		return data.reads_start_postions(min(reads_in_chain));
	}

	//TODO note duplicate in genotype_optimizer
	t_count haplo_chain_end(t_haplochain const chain) const {
		t_indices reads_in_chain(find(haplo == chain));
		t_positions end_pos = data.reads_end_postions(reads_in_chain);
		return max(end_pos);
	}

public:

	haplo_chain_optimizer(AlgorithmConfiguration const& config,
			alignment_data const& data, t_seq_bases const& ref, t_seq_bases const& alt);

	template<typename abort_checker>
	void run(abort_checker const& ac);

	template<typename abort_checker>
	void run_paired(abort_checker const& ac);

	t_haplotype get_haplotype_chains() const {
		return haplo;
	}

	t_strands get_chain_strands() const;

};

//TODO make it possible to switch between max_logsum + prior and logsum
inline t_loglike_vector haplo_chain_optimizer::compute_logsum(
		t_indices const& reads, t_position const pos) const {

	TIMER_START

	t_loglike_vector max_loglike(2);
	max_loglike(0) = -1 * std::numeric_limits<float>::infinity();
	max_loglike(1) = -1 * std::numeric_limits<float>::infinity();

	for (t_epi_base g = 0; g < 6; ++g) {

		//Ref prior
		double prior = g % 4 + 1 == ref(pos) || g % 4 + 1 == alt(pos) ? ref_match_prior_log : ref_no_match_prior_log;

		double loglike_fwd = prior;
		double loglike_rev = prior;

		for (t_index i = 0; i < reads.n_elem; ++i) {

			loglike_fwd += data.loglike_terms(reads(i), 0)(
					pos - read_start_postion(reads(i)), g);
			loglike_rev += data.loglike_terms(reads(i), 1)(
					pos - read_start_postion(reads(i)), g);
		}

		if (loglike_fwd > max_loglike(0)) {
			max_loglike(0) = loglike_fwd;
		}

		if (loglike_rev > max_loglike(1)) {
			max_loglike(1) = loglike_rev;
		}
	}

	if (!max_loglike.is_finite()) {
		throw std::runtime_error("compute_max_logsum : result not finite");
	}

	return max_loglike;

}

// DO NOT REMOVE
//inline t_loglike_vector haplo_chain_optimizer::compute_logsum(
//		t_indices const& reads, t_position const pos) const {
//
//	t_loglike_vector logsum(2, arma::fill::zeros);
//
//	for (t_epi_base g = 0; g < 6; ++g) {
//
//		double prod_fwd = 1;
//		double prod_rev = 1;
//
//		for (t_index i = 0; i < reads.n_elem; ++i) {
//
//			prod_fwd *= data.like_terms(reads(i), 0)(
//					pos - read_start_postion(reads(i)), g);
//			prod_rev *= data.like_terms(reads(i), 1)(
//					pos - read_start_postion(reads(i)), g);
//		}
//
//		logsum(0) += prod_fwd;
//		logsum(1) += prod_rev;
//
//	}
//
//	logsum = log(logsum);
//
//	if (!logsum.is_finite()) {
//		throw std::runtime_error("compute_logsum : result not finite");
//	}
//
//	return logsum;
//
//}


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

//TODO if we use logsum instead of max_logsum then we do not need ref and alt
inline haplo_chain_optimizer::haplo_chain_optimizer(
		AlgorithmConfiguration const& config, alignment_data const& data, t_seq_bases const& ref, t_seq_bases const& alt) :
				max_iterations(config.max_iterations), min_overlap_length(config.min_overlap_length),
				haplo(data.n_reads), data(data), ref(ref), alt(alt), haplochain_log_prior(config.haplochain_log_prior) {

	///////////////// Initialize:

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
		change_count = optimize_profile(ac);
		change_count += optimize_haplochains(ac);
	}

}

template<typename abort_checker>
inline void haplo_chain_optimizer::run_paired(const abort_checker& ac) {

	TIMER_START

	t_count change_count = 1;

	while(change_count != 0) {
		change_count = optimize_pair_profile(ac);
		change_count += optimize_haplochains(ac);
	}


}


//Find overlapping haplo chains
inline t_haplotype haplo_chain_optimizer::compute_feasible_haplotypes(t_position const start, t_position const end) const {

	t_index const read_number_start = min(data.read_numbers(start));
	t_index const read_number_end = max(data.read_numbers(end));

	t_haplotype feasible_haplotypes;

	for (t_index read = read_number_start; read <= read_number_end; read++) {

		//Check that read overlap haplo chain
		if (start <= read_start_postion(read)) {

			t_position overlap_length = 1 + min(end, read_end_postion(read)) - max(start, read_start_postion(read));

			//TODO this is a debug check -- should not be active in release versions
			if(overlap_length < 1) {
				//TODO remove
				//cout << overlap_length << endl;
				//cout << read_end << " : " << read_end_postion(read) << " : " << read_start << " : " << read_start_postion(read) << endl;
				throw std::runtime_error("Internal error");
			}

			if(overlap_length >= min_overlap_length) {
				feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
				feasible_haplotypes(feasible_haplotypes.n_elem - 1) = haplo(read);
			}
		}

	}

	return unique(feasible_haplotypes);
}

template<typename abort_checker>
inline t_count haplo_chain_optimizer::optimize_profile(
		const abort_checker& ac) {

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

			t_loglike_vector loglike = profile_posterior(read_number, feasible_haplotypes);

			if (!is_finite(loglike)) {
				throw std::runtime_error("optimize_profile - internal error");
			}

			t_index i = argmax(loglike);

			//check if any improvement
			if (1e-5 >= loglike(i)) { //TODO configable threshold
				continue;
			}

			//Get current haplotype chain of read
			t_haplochain old_chain = haplo(read_number);

			//Update haplotype chain of read
			haplo(read_number) = feasible_haplotypes(i);

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

			t_loglike_vector loglike ;
			loglike = profile_posterior(read_pairs, reads_pair, feasible_haplotypes);


			if (!is_finite(loglike)) {
				throw std::runtime_error("optimize_pair_profile - internal error");
			}

			t_index i = argmax(loglike);

			//check if any improvement
			if (1e-5 >= loglike(i)) { //TODO configable threshold
				continue;
			}

			//Get current haplotype chain of read pair
			t_haplochain old_chain = haplo(read_number_1);

			//Update haplotype chain of read pair
			haplo(read_number_1) = feasible_haplotypes(i);
			haplo(read_number_2) = feasible_haplotypes(i);

			//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
			update_haplotype_chain(read_pairs, haplo, old_chain);

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
inline t_count haplo_chain_optimizer::optimize_haplochains(
		const abort_checker& ac) {

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

			t_loglike_vector posterior = profile_posterior(reads_in_chain, feasible_haplotypes);

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


//TODO remove
//double haplo_chain_optimizer::compute_prior(t_haplotype const& h, t_position const start, t_position const end) const {
//
//	TIMER_START
//
//	double prior = 0;
//
//	for (t_position pos = start; pos <= end; ++pos) {
//
//			if (data.read_numbers(pos).n_elem == 0) {
//				continue;
//			}
//
//			t_haplotype unique_chains = unique(h(data.read_numbers(pos)));
//			vec nreads_chains(unique_chains.n_elem);
//			vec nreads_pos(unique_chains.n_elem);
//
//			for(u32 i = 0; i < unique_chains.n_elem; ++i) {
//				nreads_pos(i) = sum(h(data.read_numbers(pos)) == unique_chains(i));
//				nreads_chains(i) = sum(h == unique_chains(i));  // reads in chain + reads in chain over position
//			}
//
//			uvec idx = sort_index(nreads_pos, "descend");
//			nreads_pos = nreads_pos(idx);
//			nreads_chains = nreads_chains(idx);
//			idx = stable_sort_index(nreads_chains, "descend");
//
//			for(u32 i = 0; i < idx.n_elem; ++i) {
//				prior += nreads_pos(idx(i))/static_cast<double>(data.read_numbers(pos).n_elem)*compute_haplochain_prior(i+1);
//			}
//		}
//
//	return prior;
//}

double haplo_chain_optimizer::compute_chain_loglike(t_haplotype const & h,
		t_haplochain const chain) const {

	TIMER_START

	t_indices reads_in_chain = find(h == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	//find start, end pos of chain
	t_position chain_start = read_start_postion(min(reads_in_chain));
	t_position chain_end = read_end_postion(reads_in_chain(0));

	for (t_index j = 1; j < reads_in_chain.n_elem; ++j) {
		if (read_end_postion(reads_in_chain(j)) > chain_end) {
			chain_end = read_end_postion(reads_in_chain(j));
		}
	}

	t_loglike_vector logsum(2, arma::fill::zeros);

	for (t_position pos = chain_start; pos <= chain_end; ++pos) {

		t_indices reads = get_overlapping_reads_in_chain(h, pos, chain);

		//Vector of size 2 logsum of (fwd, rev)
		logsum += compute_logsum(reads, pos);

	}

	return max(logsum);
}

int haplo_chain_optimizer::chain_strand(
		t_haplochain const chain) const {

	t_indices reads_in_chain = find(haplo == chain);

	if (reads_in_chain.n_elem == 0) {
		return 0;
	}

	//find start, end pos of chain
	t_position chain_start = read_start_postion(min(reads_in_chain));
	t_position chain_end = read_end_postion(reads_in_chain(0));

	for (t_index j = 1; j < reads_in_chain.n_elem; ++j) {
		if (read_end_postion(reads_in_chain(j)) > chain_end) {
			chain_end = read_end_postion(reads_in_chain(j));
		}
	}

	t_loglike_vector logsum(2, arma::fill::zeros);

	for (t_position pos = chain_start; pos <= chain_end; ++pos) {

		t_indices reads = get_overlapping_reads_in_chain(haplo, pos, chain);

		//Vector of size 2 logsum of (fwd, rev)
		logsum += compute_logsum(reads, pos);

	}

	if(logsum(0) > logsum(1)) {
		return strand_fwd;
	}

	return strand_rev;
}

t_strands haplo_chain_optimizer::get_chain_strands() const {

	t_strands strands(haplo.n_elem);

	for(t_index i = 0; i < haplo.n_elem; ++i) {
		strands(i) = chain_strand(haplo(i));
	}

	return strands;
}


double haplo_chain_optimizer::compute_posterior(t_haplotype const& h) const {

	TIMER_START

	//Compute prior

	double posterior = 0;

	// Compute likelihood
	t_haplotype unique_chains = unique(h);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {
		posterior += compute_prior_chain(h, unique_chains(i));
		posterior += compute_chain_loglike(h, unique_chains(i));

	}

	return posterior;
}

double haplo_chain_optimizer::compute_delta_posterior(t_index const read_number, t_haplochain const new_chain) const {

	TIMER_START


	t_haplotype new_haplo(haplo); //TODO this is alot of copying
	t_haplochain current_chain = haplo(read_number);

	if(current_chain == new_chain) {
		return 0;
	}

	double delta_posterior = 0;

	//Compute posterior
	new_haplo(read_number) = new_chain;
	update_haplotype_chain(new_haplo, current_chain);

	delta_posterior += compute_prior_chain(new_haplo, new_chain);
	delta_posterior -= compute_prior_chain(haplo, new_chain);

	delta_posterior += compute_prior_chain(new_haplo, current_chain);
	delta_posterior -= compute_prior_chain(haplo, current_chain);

			//Compute delta loglike
	delta_posterior += compute_chain_loglike(new_haplo, new_chain);
	delta_posterior -= compute_chain_loglike(haplo, new_chain);

	delta_posterior += compute_chain_loglike(new_haplo, current_chain);
	delta_posterior -= compute_chain_loglike(haplo, current_chain);

	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		delta_posterior += compute_prior_chain(new_haplo, chain);
		delta_posterior += compute_chain_loglike(new_haplo, chain);
	}

	//TODO debug guards
//	if(fabs(delta_posterior - compute_posterior(new_haplo) + compute_posterior(haplo)) > 1e-5) {
//
//				cout << current_chain << " -> " << new_chain << endl;
//
//				cout << delta_posterior  << " : " << compute_posterior(new_haplo) - compute_posterior(haplo) << endl;
//				cout << compute_posterior(new_haplo)  << endl;
//				cout << compute_posterior(haplo) << endl;
//
//				for (t_index i = 0; i <= max(new_haplo); ++i) {
//					if(compute_prior_chain(new_haplo, i) != compute_prior_chain(haplo, i)) {
//						cout << "prior : " << i
//								<< " : " << compute_prior_chain(new_haplo, i)
//								<< " : " << compute_prior_chain(haplo, i) << endl;
//					}
//
//					if(compute_chain_loglike(new_haplo, i) != compute_chain_loglike(haplo, i)) {
//						cout << "loglike : " << i
//								<< " : " << compute_chain_loglike(new_haplo, i)
//								<< " : " << compute_chain_loglike(haplo, i) << endl;
//					}
//				}
//
//		throw std::runtime_error("compute_delta_posterior 1: Error");
//	}

	return delta_posterior;
}

//Merge old_chain = (reads_in_chain) and new_chain
//TODO rename - chain merge
double haplo_chain_optimizer::compute_delta_posterior(t_indices const& reads, t_haplochain const new_chain) const {

	TIMER_START

	t_haplotype new_haplo(haplo);
	t_haplochain current_chain = haplo(reads(0));

	if(current_chain == new_chain) {
		return 0;
	}

	double posterior = 0;

	//current posterior
	posterior -= compute_prior_chain(haplo, current_chain);
	posterior -= compute_prior_chain(haplo, new_chain);
	posterior -= compute_chain_loglike(haplo, current_chain);
	posterior -= compute_chain_loglike(haplo, new_chain);

	//new posterior
	new_haplo(reads).fill(new_chain);
	update_haplotype_chain(new_haplo, current_chain);

	posterior += compute_prior_chain(new_haplo, current_chain);
	posterior += compute_prior_chain(new_haplo, new_chain);
	posterior += compute_chain_loglike(new_haplo, current_chain);
	posterior += compute_chain_loglike(new_haplo, new_chain);

	//add likelihood for new chains
	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		posterior += compute_prior_chain(new_haplo, chain);
		posterior += compute_chain_loglike(new_haplo, chain);
	}

	//TODO debug guards
//	if(fabs(posterior - compute_posterior(new_haplo) + compute_posterior(haplo)) > 1e-5) {
//
//		throw std::runtime_error("compute_delta_posterior 2: Error");
//	}


	return posterior;
}

double haplo_chain_optimizer::compute_delta_posterior(
		std::vector<t_indices> read_pairs,
		t_indices const& pair,
		t_haplochain const new_chain) const {

	TIMER_START

	t_haplotype new_haplo(haplo);
	t_haplochain current_chain = haplo(pair(0));

	if(current_chain == new_chain) {
		return 0;
	}

	double posterior = 0;

	//current posterior
	posterior -= compute_prior_chain(haplo, current_chain);
	posterior -= compute_prior_chain(haplo, new_chain);
	posterior -= compute_chain_loglike(haplo, current_chain);
	posterior -= compute_chain_loglike(haplo, new_chain);

	//new posterior
	new_haplo(pair).fill(new_chain);
	update_haplotype_chain(read_pairs, new_haplo, current_chain);

	posterior += compute_prior_chain(new_haplo, current_chain);
	posterior += compute_prior_chain(new_haplo, new_chain);
	posterior += compute_chain_loglike(new_haplo, current_chain);
	posterior += compute_chain_loglike(new_haplo, new_chain);

	//add likelihood for new chains
	for(t_haplochain chain = max(haplo) + 1; chain <= max(new_haplo); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		posterior += compute_prior_chain(new_haplo, chain);
		posterior += compute_chain_loglike(new_haplo, chain);
	}

	//TODO debug guards
//	if(fabs(posterior - compute_posterior(new_haplo) + compute_posterior(haplo)) > 1) {
//
//		throw std::runtime_error("compute_delta_posterior 3: Error");
//	}


	return posterior;
}


t_loglike_vector haplo_chain_optimizer::profile_posterior(
		t_index const read_number,
		t_haplotype const& feasible_haplotypes) const {

	t_loglike_vector loglike(feasible_haplotypes.n_elem, arma::fill::zeros);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(i) = compute_delta_posterior(read_number, feasible_haplotypes(i));
	}

	return loglike;
}

t_loglike_vector haplo_chain_optimizer::profile_posterior(t_indices const& reads_in_chain, t_haplotype const& feasible_haplotypes) const {

	t_loglike_vector loglike(feasible_haplotypes.n_elem, arma::fill::zeros);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(i) = compute_delta_posterior(reads_in_chain, feasible_haplotypes(i));
	}

	return loglike;
}

t_loglike_vector haplo_chain_optimizer::profile_posterior(std::vector<t_indices> read_pairs, t_indices const& pair, t_haplotype const& feasible_haplotypes) const {

	t_loglike_vector loglike(feasible_haplotypes.n_elem, arma::fill::zeros);

	for (t_index i = 0; i < feasible_haplotypes.n_elem; i++) {
		loglike(i) = compute_delta_posterior(read_pairs, pair, feasible_haplotypes(i));
	}

	return loglike;
}


inline t_indices haplo_chain_optimizer::get_overlapping_reads_in_chain(t_haplotype const& h,
		t_position const pos, t_haplochain const chain) const {

	t_indices reads_at_pos = data.read_numbers(pos);
	arma::uvec chain_indicator = (h(reads_at_pos) == chain);

	if (sum(chain_indicator) == 0) {
		//Return empty vector
		return t_indices();
	}

	//Return read numbers of reads in chain
	return reads_at_pos(find(chain_indicator));
}

inline t_position haplo_chain_optimizer::read_start_postion(
		t_index const read_number) const {
	return data.reads_start_postions(read_number);
}

inline t_position haplo_chain_optimizer::read_end_postion(
		t_index const read_number) const {
	return data.reads_end_postions(read_number);
}

#endif /* HAPLO_CHAIN_OPTIMIZER_HPP_ */
