#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_

class haplotype {

private:

	//May only be update through move_read
	t_haplochains haplo; //vector of size n_reads
	t_strands read_strands; //vector of size n_reads
	t_positions chain_start;
	t_positions chain_end;

	//Data
	alignment_data const& data;
	std::vector<t_indices> read_pairs;
	t_positions pair_start_positions; //vector of size n_reads,
	t_positions pair_end_positions; //vector of size n_reads

	//Refs
	t_seq_bases const& ref;
	t_seq_bases const& alt;

	//Priors
	double ref_no_match_prior_log;
	double ref_match_prior_log;

	arma::Col<double> haplochain_log_prior;

	field<vec> ref_priores;


	//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
	void update_haplotype_chain(
			t_haplochains & h,
			t_haplochain const old_chain) const;

	void update_haplotype_chain(
			std::vector<t_indices> read_pairs,
			t_haplochains & h,
			t_haplochain const old_chain) const;

	void move_read(
			t_index id,
			t_strand strand,
			t_haplochain to);

	void move_read_pair(
			t_indices const& pair,
			t_strand strand,
			t_haplochain to);

	double compute_delta_posterior(
			t_index const read_number,
			t_strand strand,
			t_haplochain const new_chain) const;


	double compute_delta_posterior(
			std::vector<t_indices> read_pairs,
			t_indices const& pair,
			t_strand const strand,
			t_haplochain const new_chain) const;

	double compute_chain_loglike(
			t_haplochains const& h,
			t_strands const& strands,
			t_haplochain const chain) const;

	double compute_chain_loglike(
			t_haplochains const& h,
			t_strands const& strands,
			t_haplochain const chain,
			t_indices const& reads_overlapping_region,
			t_position region_start,
			t_position region_end) const;

	double compute_logsum(
			mat const& loglike_term,
			t_position pos) const;

	void compute_genotype(
			mat const& loglike_term,
			t_position pos,
			t_epi_base & g1,
			t_epi_base & g2) const;

	double compute_prior_chain(
			t_haplochains const& h,
			t_haplochain chain) const;

	double compute_posterior(
			t_haplochains const& h,
			t_strands const& strands) const;

	field<vec> compute_ref_priors() const;

	void init(double ref_prior);

public:

	bool use_paired_reads;

	t_count const n_elements;

	t_position const min_overlap_length;

	haplotype(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			t_position min_overlap_length,
			arma::Col<double> haplochain_log_prior,
			std::vector<t_indices> const& read_pairs);

	haplotype(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			t_position min_overlap_length,
			arma::Col<double> haplochain_log_prior);

	void move(
			t_index id,
			t_strand strand,
			t_haplochain to);

	double delta_posterior(
			t_index id,
			t_strand strand,
			t_haplochain to) const;


	t_haplochains compute_feasible_haplotypes(t_index id) const;

	void chain_clean();

	double posterior() const;

	t_genotype compute_chain_genotype(t_haplochain const chain) const;

	//getters
	field<t_genotype> genotypes() const;

	t_haplochains haplotypes() const {
		return haplo;
	}

	t_strands strands() const {
		return read_strands;
	}

	t_positions starts() const {
		return chain_start;
	}

	t_positions ends() const {
		return chain_end;
	}
};

void haplotype::init(double ref_prior) {

	DEBUG_ENTER

	chain_start.set_size(n_elements);
	chain_end.set_size(n_elements);

	if(use_paired_reads) {

		pair_start_positions.set_size(data.n_reads);
		pair_end_positions.set_size(data.n_reads);

		for(t_index i = 0; i < read_pairs.size(); ++i) {

			// assign each read pair to a chain
			haplo(read_pairs[i]).fill(i);

			// set start and end of chain
			chain_start(i) = min(data.reads_start_positions(read_pairs[i]));
			chain_end(i) = max(data.reads_end_positions(read_pairs[i]));

			// set pair start position
			pair_start_positions(read_pairs[i]).fill(chain_start(i));
			pair_end_positions(read_pairs[i]).fill(chain_end(i));

		}

	} else {

		for(t_index i = 0; i < haplo.n_elem; ++i) {

			//assign each read to a chain
			haplo(i) = i;

			// set start and end of chain
			chain_start(i) = data.reads_start_positions(i);
			chain_end(i) = data.reads_end_positions(i);

		}
	}

	//Init ref prior
	ref_no_match_prior_log = log(1-ref_prior);
	ref_match_prior_log = log(ref_prior);

	ref_priores = compute_ref_priors();

}

haplotype::haplotype(
		AlgorithmConfiguration const& config,
		const alignment_data& data,
		const t_seq_bases& ref,
		const t_seq_bases& alt,
		t_position min_overlap_length,
		arma::Col<double> haplochain_log_prior,
		const std::vector<t_indices>& read_pairs) :
				haplo(data.n_reads, fill::zeros),
				read_strands(data.n_reads, fill::zeros),
				data(data),
				read_pairs(read_pairs),
				ref(ref),
				alt(alt),
				haplochain_log_prior(config.haplochain_log_prior),
				use_paired_reads(true),
				n_elements(read_pairs.size()),
				min_overlap_length(config.min_overlap_length)
			 {

	init(config.ref_prior);

}

haplotype::haplotype(
		AlgorithmConfiguration const& config,
		const alignment_data& data,
		const t_seq_bases& ref,
		const t_seq_bases& alt,
		t_position min_overlap_length,
		arma::Col<double> haplochain_log_prior) :
				haplo(data.n_reads, fill::zeros),
				read_strands(data.n_reads, fill::zeros),
				data(data),
				read_pairs(null_read_pairs),
				ref(ref),
				alt(alt),
				haplochain_log_prior(config.haplochain_log_prior),
				use_paired_reads(false),
				n_elements(data.n_reads),
				min_overlap_length(config.min_overlap_length)
			 {
	init(config.ref_prior);
}

void haplotype::move(
		t_index id,
		t_strand strand,
		t_haplochain to) {

	if(use_paired_reads) {
		move_read_pair(read_pairs[id], strand, to);
		return;
	}

	move_read(id, strand, to);
}

double haplotype::delta_posterior(
		t_index id,
		t_strand strand,
		t_haplochain to) const {

	if(use_paired_reads) {
		return compute_delta_posterior(read_pairs, read_pairs[id], strand, to);
	}

	return compute_delta_posterior(id, strand, to);

}

void haplotype::move_read(
		t_index id,
		t_strand strand,
		t_haplochain to) {

	TIMER_START
	DEBUG_ENTER

	t_haplochain from = haplo(id);

	if(from == to && strand == read_strands(id)) {
		return;
	}

	//Update haplo and strand
	t_haplochain max_chain_old = max(haplo);

	haplo(id) = to;
	read_strands(id) = strand;

	update_haplotype_chain(haplo, from);

	//update chain start and end
	if(to <= max_chain_old) {
		chain_start(to) = min(chain_start(to), data.reads_start_positions(id));
		chain_end(to) = max(chain_end(to), data.reads_end_positions(id));
	}

	t_indices reads_in_chain(find(haplo == from)); //TODO this could properly be done more efficient
	if(reads_in_chain.n_elem > 0) {
		chain_start(from) = min(data.reads_start_positions(reads_in_chain));
    	chain_end(from) =  max(data.reads_end_positions(reads_in_chain));
	}

	//Add new chains if any
	chain_start.resize(max(haplo)+1);
	chain_end.resize(max(haplo)+1);

	for(t_haplochain chain = max_chain_old + 1; chain <= max(haplo); ++chain) {

		if(chain == from) {
			continue;
		}

		t_indices reads(find(haplo == chain)); //TODO this could properly be done more efficient
		chain_start(chain) = min(data.reads_start_positions(reads));
	    chain_end(chain) =  max(data.reads_end_positions(reads));

	}
}


void haplotype::move_read_pair(
		t_indices const& pair,
		t_strand strand,
		t_haplochain to) {

	TIMER_START
	DEBUG_ENTER

	t_haplochain from = haplo(pair(0));

	if(from == to && strand == read_strands(pair(0))) {
		return;
	}

	//Update haplo and strand
	t_haplochain max_chain_old = max(haplo);

	haplo(pair).fill(to);
	read_strands(pair).fill(strand);

	update_haplotype_chain(read_pairs, haplo, from);

	//update chain start and end
	if(to <= max_chain_old) {
		chain_start(to) = min(chain_start(to), min(data.reads_start_positions(pair)));
		chain_end(to) = max(chain_end(to), max(data.reads_end_positions(pair)));
	}

	t_indices reads_in_chain(find(haplo == from)); //TODO this could properly be done more efficient
	if(reads_in_chain.n_elem > 0) {
		chain_start(from) = min(data.reads_start_positions(reads_in_chain));
		chain_end(from) =  max(data.reads_end_positions(reads_in_chain));
	}

	//Add new chains if any
	chain_start.resize(max(haplo)+1);
	chain_end.resize(max(haplo)+1);

	for(t_haplochain chain = max_chain_old + 1; chain <= max(haplo); ++chain) {

		if(chain == from) {
			continue;
		}

		t_indices reads(find(haplo == chain)); //TODO this could properly be done more efficient
		chain_start(chain) = min(data.reads_start_positions(reads));
	    chain_end(chain) =  max(data.reads_end_positions(reads));

	}
}


void haplotype::update_haplotype_chain(
		t_haplochains & haplotype,
		t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_indices reads_in_chain = find(haplotype == old_chain);

	t_haplochain chain = old_chain;
	t_position chain_end = data.reads_end_positions(reads_in_chain(0));

	for (t_index i = 1; i < reads_in_chain.n_elem; ++i) {

		t_index read = reads_in_chain(i);

		t_position overlap_length = chain_end - data.reads_start_positions(read) + 1;
		if(overlap_length < min_overlap_length) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;

			haplotype(read) = chain;
			chain_end = data.reads_end_positions(read);
		}

		else {
			haplotype(read) = chain;
			chain_end = max(chain_end, data.reads_end_positions(read));
		}
	}
}

void haplotype::update_haplotype_chain(
		std::vector<t_indices> read_pairs,
		t_haplochains & haplotype,
		t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_haplochain chain = old_chain;
	t_position chain_end = max(data.reads_end_positions(read_pairs[0](0)), data.reads_end_positions(read_pairs[0](1)));

	for (t_index i = 1; i < read_pairs.size(); ++i) {

		//check if read pair is in chain
		if(haplotype(read_pairs[i](0)) != old_chain) {
			continue;
		}

		t_position start = min(data.reads_start_positions(read_pairs[i](0)), data.reads_start_positions(read_pairs[i](1)));

		if (start > chain_end - min_overlap_length + 1) {
			//reads not overlapping => new chain
			chain = max(haplotype) + 1;

			haplotype(read_pairs[i](0)) = chain;
			haplotype(read_pairs[i](1)) = chain;

			chain_end = max(data.reads_end_positions(read_pairs[i](0)), data.reads_end_positions(read_pairs[i](1)));

		}

		else {
			haplotype(read_pairs[i](0)) = chain;
			haplotype(read_pairs[i](1)) = chain;

			t_position read_pair_end = max(data.reads_end_positions(read_pairs[i](0)), data.reads_end_positions(read_pairs[i](1)));
			chain_end = max(chain_end, read_pair_end);
		}
	}
}

void haplotype::chain_clean() {

	t_haplochains unique_chains = sort(unique(haplo));

	t_haplochains new_haplo(haplo);
	t_positions new_chain_start(unique_chains.n_elem);
	t_positions new_chain_end(unique_chains.n_elem);

	for (t_index i = 0; i < unique_chains.n_elem; ++i) {
		t_indices reads_in_chain(find(haplo == unique_chains(i)));
		new_haplo(reads_in_chain).fill(i);
		new_chain_start(i) = chain_start(unique_chains(i));
		new_chain_end(i) = chain_end(unique_chains(i));
	}

	haplo = new_haplo;
	chain_start = new_chain_start;
	chain_end = new_chain_end;
}



double haplotype::compute_delta_posterior(
		t_index const read_number,
		t_strand strand,
		t_haplochain const new_chain) const {

	TIMER_START
	DEBUG_ENTER

	t_haplochains new_haplo(haplo); //TODO this is alot of copying
	t_strands new_strands(read_strands); //TODO this is alot of copying

	t_haplochain current_chain = haplo(read_number);

	if(current_chain == new_chain && strand == read_strands(read_number)) {
		return 0;
	}

	t_position region_start = data.reads_start_positions(read_number)-1;
	t_position region_end = data.reads_end_positions(read_number)+1;

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

double haplotype::compute_delta_posterior(
		std::vector<t_indices> read_pairs,
		t_indices const& pair,
		t_strand const strand,
		t_haplochain const new_chain) const {

	TIMER_START

	t_haplochains new_haplo(haplo);
	t_strands new_strands(read_strands); //TODO this is alot of copying

	t_haplochain current_chain = haplo(pair(0));

	if(current_chain == new_chain && strand == read_strands(pair(0))) {
		return 0;
	}

	t_position region_start = min(data.reads_start_positions(pair))-1;
	t_position region_end = max(data.reads_end_positions(pair))+1;

	t_indices reads_in_region = find(pair_start_positions <= region_end && pair_end_positions >= region_start);


	//new posterior
	new_haplo(pair).fill(new_chain);
	new_strands(pair).fill(strand);
	update_haplotype_chain(read_pairs, new_haplo, current_chain);

	double delta_posterior = 0;

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
//		throw std::runtime_error("compute_delta_posterior 3: Error");
//	}


	return delta_posterior;
}


double haplotype::compute_chain_loglike(
		t_haplochains const& h,
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
		t_position const read_start = data.reads_start_positions(*r);
		t_position const read_end = data.reads_end_positions(*r);

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


double haplotype::compute_chain_loglike(
		t_haplochains const& h,
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
		t_position const read_start = data.reads_start_positions(*r);
		t_position const read_end = data.reads_end_positions(*r);

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


double haplotype::compute_logsum(
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


field<vec> haplotype::compute_ref_priors() const {

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

t_genotype haplotype::compute_chain_genotype(t_haplochain const chain) const {

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
		t_position const read_start = data.reads_start_positions(*r);
		t_position const read_end = data.reads_end_positions(*r);

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

void haplotype::compute_genotype(
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


field<t_genotype> haplotype::genotypes() const {

	t_haplochains unique_chains = sort(unique(haplo));

	field<t_genotype> genotypes(unique_chains.n_elem);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			genotypes(i) = compute_chain_genotype(unique_chains(i));
		}

	return(genotypes);
}

t_haplochains haplotype::compute_feasible_haplotypes(t_index id) const {

	DEBUG_ENTER
	TIMER_START

	t_haplochains feasible_haplotypes;

	t_position start;
	t_position end;

	if(use_paired_reads) {
		t_indices pair = read_pairs[id];
		start = min(data.reads_start_positions(pair));
		end = max(data.reads_end_positions(pair));
	}

	else {
		start = data.reads_start_positions(id);
		end = data.reads_end_positions(id);
	}

	//Compute for each chain the number of bases overlapping
	t_haplochains unique_chains = unique(haplo);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {

		t_position overlap_length = pos(end-start + 1 - pos(chain_start(unique_chains(i)) - start) - pos(end-chain_end(unique_chains(i))));

		if(overlap_length >= min_overlap_length) {
			feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
			feasible_haplotypes(feasible_haplotypes.n_elem - 1) = unique_chains(i);
		}

	}

	//Add free haplo chain , i.e. a haplo chain with no reads
	feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
	feasible_haplotypes(feasible_haplotypes.n_elem - 1) = max(haplo) + 1;

	return feasible_haplotypes;

}

double haplotype::compute_prior_chain(
		t_haplochains const& h,
		t_haplochain chain) const {

	TIMER_START

	t_count reads_in_chain = sum(h == chain);

	if(reads_in_chain == 0) {
		return 0;
	}

	return haplochain_log_prior(reads_in_chain-1);

}

double haplotype::compute_posterior(
		t_haplochains const& h,
		t_strands const& strands) const {

	TIMER_START

	//Compute prior

	double posterior = 0;

	// Compute likelihood
	t_haplochains unique_chains = unique(h);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {
		posterior += compute_prior_chain(h, unique_chains(i));
		posterior += compute_chain_loglike(h, strands, unique_chains(i));

	}

	return posterior;
}

double haplotype::posterior() const {
	return compute_posterior(haplo, read_strands);
}

#endif /* HAPLOTYPE_HPP_ */
