#ifndef HAPLOTYPE_HPP_
#define HAPLOTYPE_HPP_

class haplotype_state {

public:

	t_haplochains haplo; //vector of size n_reads
	t_strands read_strands; //vector of size n_reads
	t_positions chain_start; //vector of size n chais
	t_positions chain_end; //vector of size n chais
	double prior;

	haplotype_state(int n_reads) :
		haplo(n_reads, fill::zeros),
		read_strands(n_reads, fill::zeros) {
	}

};

class haplotype {

private:

	//May only be update through move_read
	haplotype_state state;

	//Data
	alignment_data const& data;

	//Blocks
	std::vector<t_indices> const read_blocks;
	t_positions const block_start_positions; //vector of size n_reads,
	t_positions const block_end_positions; //vector of size n_reads

	//Refs
	t_seq_bases const& ref;
	t_seq_bases const& alt;

	//Priors
	double const ref_no_match_prior_log;
	double const ref_match_prior_log;

	field<vec> const ref_priores;

	bool is_read_feasible(
			t_position read_start,
			t_position read_end,
			t_position chain_start,
			t_position chain_end) const;

	bool is_block_feasible(
			t_indices const& block,
			t_position chain_start,
			t_position chain_end) const;

	//Update old haplotype chain to ensure that it is still feasible (after removal of read_number)
	void update_haplotype_chain(
			t_haplochains & h,
			t_haplochain const old_chain) const;

	void update_haplotype_chain(
			std::vector<t_indices> const& read_blocks,
			t_haplochains & h,
			t_haplochain const old_chain) const;

	//compute blocks of reads which are feasible
	std::vector<t_indices> feasible_blocks(t_indices const& reads) const;

	std::vector<t_indices> get_feasible_blocks(t_indices const& reads) const;

	std::vector<t_indices> get_feasible_blocks(std::vector<t_indices> const& rblocks) const;


	void move_read(
			t_haplochains & h,
			t_strands & s,
			t_positions & c_start,
			t_positions & c_end,
			double & prior,
			t_index id,
			t_strand strand,
			t_haplochain to) const;

	void move_read_block(
			t_haplochains & h,
			t_strands & s,
			t_positions & c_start,
			t_positions & c_end,
			double & prior,
			t_indices const& block,
			t_strand strand,
			t_haplochain to) const;

	double compute_delta_posterior(
			t_index const read_number,
			t_strand strand,
			t_haplochain const new_chain) const;


	double compute_delta_posterior(
			std::vector<t_indices> read_blocks,
			t_indices const& block,
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
			t_positions const& c_start,
			t_positions const& c_end,
			t_haplochain chain) const;

	double compute_posterior(
			t_haplochains const& h,
			t_positions const& c_start,
			t_positions const& c_end,
			t_strands const& strands) const;

	double compute_prior(
			t_haplochains const& h,
			t_positions const& c_start,
			t_positions const& c_end,
			t_strands const& strands) const;

	double compute_delta(
			t_haplochains const& h,
			t_positions const& c_start,
			t_positions const& c_end,
			double prior,
			t_strands const& s,
			t_haplochains const& new_h,
			t_positions const& new_c_start,
			t_positions const& new_c_end,
			double new_prior,
			t_strands const& new_s,
			t_indices const& reads_in_region,
			t_position region_start,
			t_position region_end,
			t_haplochain current_chain,
			t_haplochain new_chain) const;

	field<vec> compute_ref_priors() const;

	void init(double ref_prior);

	bool is_HCGD(t_position pos) const {
		return (pos > 0 && ref(pos-1) != 2 && ref(pos) == 1 && ref(pos+1) == 2 && ref(pos+2) != 1);
	}

	bool is_DGCH(t_position pos) const {
		return (pos > 0 && ref(pos-1) != 1 && ref(pos) == 2 && ref(pos+1) == 1 && ref(pos+2) != 2);
	}

	t_count count_CpG(t_position start, t_position end) const {

	//	cout << ref.n_elem << " : " << start << " : " << end << endl;

		return sum(ref.subvec(start, end) == 1 && ref.subvec(start+1, end+1) == 2);
	}

	t_count count_HCGD(t_position start, t_position end) const {

		if(start == 0) {
			start = 1;
		}

		return sum(ref.subvec(start-1, end-1) != 2
				&& ref.subvec(start, end) == 1
				&& ref.subvec(start+1, end+1) == 2
				&& ref.subvec(start+2, end+2) != 1);
	}

	t_count count_DGCH(t_position start, t_position end) const {

		if(start == 0) {
			start = 1;
		}

		return sum(ref.subvec(start-1, end-1) != 1
				&& ref.subvec(start, end) == 2
				&& ref.subvec(start+1, end+1) == 1
				&& ref.subvec(start+2, end+2) != 2);
	}

public:

	bool const use_read_blocks;
	bool const split_mode;

	t_count const n_elements;

	t_position const min_overlap_length;
	t_count const min_CG_count;
	t_count const min_HCGD_count;
	t_count const min_DGCH_count;
	t_position const margin;

	double const prior_adjust;

	haplotype(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			t_position min_overlap_length,
			std::vector<t_indices> const& read_blocks);

	haplotype(
			AlgorithmConfiguration const& config,
			alignment_data const& data,
			t_seq_bases const& ref,
			t_seq_bases const& alt,
			t_position min_overlap_length);

	void move(
			t_index id,
			t_strand strand,
			t_haplochain to);

	double delta_posterior(
			t_index id,
			t_strand strand,
			t_haplochain to) const;

	void set_blocks(haplotype const& h);

	void set_blocks(std::vector<t_indices> const& blocks);

	void set_blocks();

	void reset_haplotype();

	void set_haplochain_prior(vec const& log_prior);

	void set_min(t_position overlap,
			t_count CG_count,
			t_count HCGD_count,
			t_count DGCH_count);

	t_haplochains compute_feasible_haplotypes(t_index id) const;

	void chain_clean();

	double posterior() const;

	t_genotype compute_chain_genotype(t_haplochain const chain) const;

	mat compute_chain_loglike(t_haplochain const chain) const;

	//getters
	field<t_genotype> genotypes() const;

	field<mat> loglikes() const;

	t_haplochains haplotypes() const {
		return state.haplo;
	}

	t_strand strand(t_index id) const {

		if(use_read_blocks) {
			return state.read_strands(read_blocks[id](0));
		}

		return state.read_strands(id);
	}

	t_strands strands() const {
		return state.read_strands;
	}

	t_positions starts() const {
		return state.chain_start;
	}

	t_positions ends() const {
		return state.chain_end;
	}

	haplotype_state get_state() const {
		return state;
	}

	void set_state(haplotype_state const& s) {
		state = s;
	}
};

std::vector<t_indices> haplotype::feasible_blocks(t_indices const& reads) const {

	if(use_read_blocks) {

		std::vector<t_indices> blocks;

		for (t_index i = 0; i < read_blocks.size(); ++i) {
			for(t_index j = 0; j < read_blocks[i].n_elem; ++j) {
				if(any(reads == read_blocks[i](j))) {
					blocks.push_back(read_blocks[i]);
					break;
				}
			}
		}

		return get_feasible_blocks(blocks);
	}

	return get_feasible_blocks(reads);
}


std::vector<t_indices> haplotype::get_feasible_blocks(t_indices const& reads) const {

	std::vector<t_indices> blocks;

	t_indices block(1);
	block(0) = reads(0);

	t_position chain_end = data.reads_end_positions(reads(0));

	for (t_index i = 1; i < reads.n_elem; ++i) {

		t_index read = reads(i);

		if(! is_read_feasible(data.reads_start_positions(read),
				data.reads_end_positions(read),
				data.reads_start_positions(read),
				chain_end)) {

			//not feasible => new chain
			blocks.push_back(block);
			block.set_size(1);
			block(0) = read;

			chain_end = data.reads_end_positions(read);
		}

		else {
			block.resize(block.n_elem + 1);
			block(block.n_elem - 1) = read;

			chain_end = max(chain_end, data.reads_end_positions(read));
		}
	}

	if(block.n_elem > 0) {
		blocks.push_back(block);
	}

	return blocks;
}

std::vector<t_indices> haplotype::get_feasible_blocks(std::vector<t_indices> const& rblocks) const {

	std::vector<t_indices> blocks;

	t_indices block(rblocks[0]);

	t_position chain_end = block_end_positions(rblocks[0](0));

	for (t_index i = 1; i < rblocks.size(); ++i) {

		t_position start = block_start_positions(rblocks[i](0));
		t_position end = block_end_positions(rblocks[i](0));

		if(! is_block_feasible(rblocks[i], start, chain_end)) {

			//not feasible => new chain
			blocks.push_back(block);
			block.reset();
			block = rblocks[i];

			chain_end = end;
		}

		else {
			block = join_cols(block, rblocks[i]);

			chain_end = max(chain_end, end);
		}
	}

	if(block.n_elem > 0) {
		blocks.push_back(block);
	}

	return blocks;
}

void haplotype::set_blocks(haplotype const& h) {

	DEBUG_ENTER

	std::vector<t_indices> blocks;
	t_haplochains const hc = h.haplotypes();
	t_strands const rs = h.strands();

	for(t_index i = 0; i <= max(hc); ++i) {

		t_indices block_fwd = find(hc == i && rs == strand_fwd);
		if(block_fwd.n_elem != 0) {

			std::vector<t_indices> tmp = feasible_blocks(block_fwd);
			blocks.insert(blocks.end(), tmp.begin(), tmp.end());

			//blocks.push_back(block_fwd);
		}

		t_indices block_rev = find(hc == i && rs == strand_rev);
		if(block_rev.n_elem != 0) {


			std::vector<t_indices> tmp = feasible_blocks(block_rev);
			blocks.insert(blocks.end(), tmp.begin(), tmp.end());

			//blocks.push_back(block_rev);
		}
	}

	set_blocks(blocks);
}

void haplotype::set_blocks(std::vector<t_indices> const& blocks) {

	//Test validity of blocks'
	//TODO debug guard
	t_positions read_present(data.n_reads, arma::fill::zeros);

	for(t_index i = 0; i < blocks.size(); ++i) {
		if(sum(read_present(blocks[i])) > 0) {
			throw std::runtime_error("set_blocks - blocks not unique error");
		}
		read_present(blocks[i]).ones();
	}

	if(any(read_present == 0)) {
		std::cout << trans(read_present) << endl;
		throw std::runtime_error("set_blocks - blocks not compete error");
	}

	//Get block positions
	t_positions block_start(blocks.size());

	t_positions block_start_positions_tmp(data.n_reads);
	t_positions block_end_positions_tmp(data.n_reads);

	for(t_index i = 0; i < blocks.size(); ++i) {

				// get start and end of block
				t_position start = min(data.reads_start_positions(blocks[i]));
				t_position end = max(data.reads_end_positions(blocks[i]));

				// set block start position
				block_start(i) = start;

				// set block positions
				block_start_positions_tmp(blocks[i]).fill(start);
				block_end_positions_tmp(blocks[i]).fill(end);

			}

	//Reorder blocks such that they are ordered with increasing start position
	uvec idx = sort_index(block_start);
	std::vector<t_indices> blocks_ordered;
	for(t_index i = 0; i < idx.n_elem; ++i) {

		//order blocks
		blocks_ordered.push_back(blocks[idx(i)]);
	}

	const_cast<bool&>(use_read_blocks) = true;
	const_cast<t_count&>(n_elements) = blocks_ordered.size();
	const_cast<std::vector<t_indices> &>(read_blocks) = blocks_ordered;

	const_cast<t_positions&>(block_start_positions) = block_start_positions_tmp;
	const_cast<t_positions&>(block_end_positions) = block_end_positions_tmp;

}

void haplotype::set_blocks() {
	const_cast<bool&>(use_read_blocks) = false;
	const_cast<t_count&>(n_elements) = data.n_reads;
}


void haplotype::set_min(t_position overlap,
		t_count CG_count,
		t_count HCGD_count,
		t_count DGCH_count) {

	const_cast<t_position&>(min_overlap_length) = overlap;
	const_cast<t_count&>(min_CG_count) = CG_count;
	const_cast<t_count&>(min_HCGD_count) = HCGD_count;
	const_cast<t_count&>(min_DGCH_count) = DGCH_count;

}

void haplotype::init(double ref_prior) {

	DEBUG_ENTER

	state.chain_start.set_size(n_elements);
	state.chain_end.set_size(n_elements);

	if(use_read_blocks) {

		t_positions block_start_positions_tmp(data.n_reads);
		t_positions block_end_positions_tmp(data.n_reads);

		for(t_index i = 0; i < read_blocks.size(); ++i) {

			// assign each read block to a chain
			state.haplo(read_blocks[i]).fill(i);

			// set start and end of chain
			state.chain_start(i) = min(data.reads_start_positions(read_blocks[i]));
			state.chain_end(i) = max(data.reads_end_positions(read_blocks[i]));

			// set block start position
			block_start_positions_tmp(read_blocks[i]).fill(state.chain_start(i));
			block_end_positions_tmp(read_blocks[i]).fill(state.chain_end(i));

		}

		 const_cast<t_positions&>(block_start_positions) = block_start_positions_tmp;
		 const_cast<t_positions&>(block_end_positions) = block_end_positions_tmp;

	} else {

		for(t_index i = 0; i < state.haplo.n_elem; ++i) {

			//assign each read to a chain
			state.haplo(i) = i;

			// set start and end of chain
			state.chain_start(i) = data.reads_start_positions(i);
			state.chain_end(i) = data.reads_end_positions(i);

		}
	}

	//Init structural prior
	state.prior = 0;
	for (t_index i = 0; i < state.haplo.n_elem; ++i) {
		//Compute prior
		state.prior += compute_prior_chain(state.haplo, state.chain_start, state.chain_end, i);
	}

	//Init ref prior
	const_cast<double&>(ref_no_match_prior_log) = log(1-ref_prior);
	const_cast<double&>(ref_match_prior_log) = log(ref_prior);

	const_cast<field<vec>&>(ref_priores) = compute_ref_priors();

}

void haplotype::reset_haplotype() {

	const_cast<t_count&>(n_elements) = read_blocks.size();

	state.chain_start.set_size(n_elements);
	state.chain_end.set_size(n_elements);

	t_positions block_start_positions_tmp(data.n_reads);
	t_positions block_end_positions_tmp(data.n_reads);

	for(t_index i = 0; i < read_blocks.size(); ++i) {

		// assign each read block to a chain
		state.haplo(read_blocks[i]).fill(i);

		// set start and end of chain
		state.chain_start(i) = min(data.reads_start_positions(read_blocks[i]));
		state.chain_end(i) = max(data.reads_end_positions(read_blocks[i]));

		// set block start position
		block_start_positions_tmp(read_blocks[i]).fill(state.chain_start(i));
		block_end_positions_tmp(read_blocks[i]).fill(state.chain_end(i));

	}

	const_cast<t_positions&>(block_start_positions) = block_start_positions_tmp;
	const_cast<t_positions&>(block_end_positions) = block_end_positions_tmp;

	//Init structural prior
	state.prior = 0;
	for (t_index i = 0; i < state.haplo.n_elem; ++i) {
		//Compute prior
		state.prior += compute_prior_chain(state.haplo, state.chain_start, state.chain_end, i);
	}

}

haplotype::haplotype(
		AlgorithmConfiguration const& config,
		const alignment_data& data,
		const t_seq_bases& ref,
		const t_seq_bases& alt,
		t_position min_overlap_length,
		const std::vector<t_indices>& read_blocks) :
				state(data.n_reads),
				data(data),
				read_blocks(read_blocks),
				ref(ref),
				alt(alt),
				ref_no_match_prior_log(),
				ref_match_prior_log(),
				ref_priores(),
				use_read_blocks(true),
				split_mode(config.split_mode),
				n_elements(read_blocks.size()),
				min_overlap_length(config.min_overlap_length),
				min_CG_count(config.min_CG_count),
				min_HCGD_count(config.min_HCGD_count),
				min_DGCH_count(config.min_DGCH_count),
				margin(config.margin),
				prior_adjust(config.structual_prior_scale)
			 {

	init(config.ref_prior);

}

haplotype::haplotype(
		AlgorithmConfiguration const& config,
		const alignment_data& data,
		const t_seq_bases& ref,
		const t_seq_bases& alt,
		t_position min_overlap_length) :
				state(data.n_reads),
				data(data),
				read_blocks(null_read_blocks),
				ref(ref),
				alt(alt),
				ref_no_match_prior_log(),
				ref_match_prior_log(),
				ref_priores(),
				use_read_blocks(false),
				split_mode(config.split_mode),
				n_elements(data.n_reads),
				min_overlap_length(config.min_overlap_length),
				min_CG_count(config.min_CG_count),
				min_HCGD_count(config.min_HCGD_count),
				min_DGCH_count(config.min_DGCH_count),
				margin(config.margin),
				prior_adjust(config.structual_prior_scale)
			 {

	init(config.ref_prior);
}

void haplotype::move(
		t_index id,
		t_strand strand,
		t_haplochain to) {

	if(use_read_blocks) {
		move_read_block(state.haplo, state.read_strands, state.chain_start, state.chain_end, state.prior, read_blocks[id], strand, to);
		return;
	}

	move_read(state.haplo, state.read_strands, state.chain_start, state.chain_end, state.prior, id, strand, to);
}

double haplotype::delta_posterior(
		t_index id,
		t_strand strand,
		t_haplochain to) const {

	if(use_read_blocks) {
		return compute_delta_posterior(read_blocks, read_blocks[id], strand, to);
	}

	return compute_delta_posterior(id, strand, to);

}

void haplotype::move_read(
		t_haplochains & h,
		t_strands & s,
		t_positions & c_start,
		t_positions & c_end,
		double & prior,
		t_index id,
		t_strand strand,
		t_haplochain to) const {

	TIMER_START
	DEBUG_ENTER

	t_haplochain from = h(id);

	if(from == to && strand == s(id)) {
		return;
	}

	//update prior
	prior -= compute_prior_chain(h, c_start, c_end, from);
	prior -= compute_prior_chain(h, c_start, c_end, to);

	//Update haplo and strand
	t_haplochain max_chain_old = max(h);

	h(id) = to;
	s(id) = strand;

	update_haplotype_chain(h, from);

	//update chain start and end
	if(to <= max_chain_old) {
		c_start(to) = min(c_start(to), data.reads_start_positions(id));
		c_end(to) = max(c_end(to), data.reads_end_positions(id));
	}

	t_indices reads_in_chain(find(h == from)); //TODO this could properly be done more efficient
	if(reads_in_chain.n_elem > 0) {
		c_start(from) = min(data.reads_start_positions(reads_in_chain));
		c_end(from) =  max(data.reads_end_positions(reads_in_chain));
	}

	//Add new chains if any
	c_start.resize(max(h)+1);
	c_end.resize(max(h)+1);

	for(t_haplochain chain = max_chain_old + 1; chain <= max(h); ++chain) {

		if(chain == from) {
			continue;
		}

		t_indices reads(find(h == chain)); //TODO this could properly be done more efficient
		c_start(chain) = min(data.reads_start_positions(reads));
		c_end(chain) =  max(data.reads_end_positions(reads));
	}

	//update prior
	for(t_haplochain chain = max_chain_old + 1; chain <= max(h); ++chain) {

		if(chain == from) {
			continue;
		}

		prior += compute_prior_chain(h, c_start, c_end, chain);
	}

	prior += compute_prior_chain(h, c_start, c_end, from);

	if(to <= max_chain_old) {
		prior += compute_prior_chain(h, c_start, c_end, to);
	}
}


void haplotype::move_read_block(
		t_haplochains & h,
		t_strands & s,
		t_positions & c_start,
		t_positions & c_end,
		double & prior,
		t_indices const& block,
		t_strand strand,
		t_haplochain to) const {

	TIMER_START
	DEBUG_ENTER

	t_haplochain from = h(block(0));

	if(from == to && strand == s(block(0))) {
		return;
	}

	//update prior
	prior -= compute_prior_chain(h, c_start, c_end, from);
	prior -= compute_prior_chain(h, c_start, c_end, to);

	//Update haplo and strand
	t_haplochain max_chain_old = max(h);

	h(block).fill(to);
	s(block).fill(strand);

	update_haplotype_chain(read_blocks, h, from);

	//update chain start and end
	if(to <= max_chain_old) {
		c_start(to) = min(c_start(to), min(data.reads_start_positions(block)));
		c_end(to) = max(c_end(to), max(data.reads_end_positions(block)));
	}

	t_indices reads_in_chain(find(h == from)); //TODO this could properly be done more efficient
	if(reads_in_chain.n_elem > 0) {
		c_start(from) = min(data.reads_start_positions(reads_in_chain));
		c_end(from) =  max(data.reads_end_positions(reads_in_chain));
	}

	//Add new chains if any
	c_start.resize(max(h)+1);
	c_end.resize(max(h)+1);

	for(t_haplochain chain = max_chain_old + 1; chain <= max(h); ++chain) {

		if(chain == from) {
			continue;
		}

		t_indices reads(find(h == chain)); //TODO this could properly be done more efficient
		c_start(chain) = min(data.reads_start_positions(reads));
		c_end(chain) =  max(data.reads_end_positions(reads));

	}

	//update prior
	for(t_haplochain chain = max_chain_old + 1; chain <= max(h); ++chain) {

		if(chain == from) {
			continue;
		}

		prior += compute_prior_chain(h, c_start, c_end, chain);
	}

	prior += compute_prior_chain(h, c_start, c_end, from);

	if(to <= max_chain_old) {
		prior += compute_prior_chain(h, c_start, c_end, to);
	}

//#ifdef CHECK_POSTERIOR
//	double p = compute_prior(h, c_start, c_end, s);
//	if(fabs(p - prior) > 1e-5) {
//		throw std::runtime_error("move_read_block Error");
//	}
//#endif

}

bool haplotype::is_read_feasible(
		t_position read_start,
		t_position read_end,
		t_position chain_start,
		t_position chain_end) const {

	t_position overlap_start = max(read_start, chain_start);
	t_position overlap_end = min(read_end, chain_end);

	if(overlap_end - overlap_start + 1 < min_overlap_length) {
		return false;
	}

	//TODO configable margin
	overlap_start = overlap_start + margin;
	overlap_end = overlap_end - margin;

	if(overlap_start < 0 || overlap_end <= overlap_start) {
		return false;
	}

//	cout << overlap_end << " : " << overlap_start << endl;
//	cout << overlap_end + 1 << " >= " << min_overlap_length + overlap_start<< endl;
//
	if(min_CG_count > 0 && count_CpG(overlap_start, overlap_end) >= min_CG_count) {
		return true;
	}

	if(min_HCGD_count > 0 && count_HCGD(overlap_start, overlap_end) >= min_HCGD_count) {
		return true;
	}

	if(min_DGCH_count > 0 && count_DGCH(overlap_start, overlap_end) >= min_DGCH_count) {
		return true;
	}

	return false;

}

bool haplotype::is_block_feasible(
		t_indices const& block,
		t_position chain_start,
		t_position chain_end) const {

	t_position block_start = block_start_positions(block(0));
	t_position block_end = block_end_positions(block(0));

	t_position block_overlap_start = max(block_start, chain_start);
	t_position block_overlap_end = min(block_end, chain_end);

	if(block_overlap_end - block_overlap_start + 1 < min_overlap_length) {
		return false;
	}

	//TODO configable margin
	block_overlap_start = block_overlap_start + margin;
	block_overlap_end = block_overlap_end - margin;

	if(block_overlap_start < 0 || block_overlap_end <= block_overlap_start) {
		return false;
	}

	//	cout << overlap_end << " : " << overlap_start << endl;
	//	cout << overlap_end + 1 << " >= " << min_overlap_length + overlap_start<< endl;

	t_count c_GC = 0;
	t_count c_HCGD = 0;
	t_count c_DGCH = 0;

	t_position end_old = 0;

	for(t_count i = 0; i < block.n_elem; ++i) {

		t_position start = data.reads_start_positions(block(i));
		t_position end = data.reads_end_positions(block(i));

		if(end <= end_old) {
			continue;
		}

		start = max(start, end_old);
		end_old = end;

		t_position overlap_start = max(start, chain_start);
		t_position overlap_end = min(end, chain_end);

		overlap_start = overlap_start + margin;
		overlap_end = overlap_end - margin;

		if(overlap_start >= 0 && overlap_end > overlap_start) {

			if(min_CG_count > 0) {
				c_GC += count_CpG(overlap_start, overlap_end);

				if(c_GC >= min_CG_count) {
					return true;
				}
			}

			if(min_HCGD_count > 0) {
				c_HCGD += count_HCGD(overlap_start, overlap_end);

				if(c_HCGD >= min_HCGD_count) {
					return true;
				}
			}

			if(min_DGCH_count > 0) {
				c_DGCH += count_DGCH(overlap_start, overlap_end);

				if(c_DGCH >= min_DGCH_count) {
					return true;
				}
			}
		}
	}

	return false;

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

		if(! is_read_feasible(data.reads_start_positions(read),
				data.reads_end_positions(read),
				data.reads_start_positions(read),
				chain_end)) {

			//not feasible => new chain
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
		std::vector<t_indices> const& read_blocks,
		t_haplochains & haplotype,
		t_haplochain const old_chain) const {

	if (sum(haplotype == old_chain) == 0) {
		//empty chain
		return;
	}

	t_haplochain chain = old_chain;
	t_position chain_end;
	bool first_read_block = true;

	for (t_index i = 1; i < read_blocks.size(); ++i) {

		//check if read block is in chain
		if(haplotype(read_blocks[i](0)) != old_chain) {
			continue;
		}

		if(first_read_block) {
			chain_end = block_end_positions(read_blocks[i](0));// max(data.reads_end_positions(read_blocks[i]));
			first_read_block = false;
			continue;
		}

		t_position start = block_start_positions(read_blocks[i](0));//min(data.reads_start_positions(read_blocks[i]));
		t_position end = block_end_positions(read_blocks[i](0));

		if(! is_block_feasible(read_blocks[i], start, chain_end)) {
			//not feasible=> new chain
			chain = max(haplotype) + 1;

			haplotype(read_blocks[i]).fill(chain);

			chain_end = max(data.reads_end_positions(read_blocks[i]));

		}

		else {
			haplotype(read_blocks[i]).fill(chain);
			chain_end = max(chain_end, end);
		}
	}
}

void haplotype::chain_clean() {

	t_haplochains unique_chains = sort(unique(state.haplo));

	t_haplochains new_haplo(state.haplo);
	t_positions new_chain_start(unique_chains.n_elem);
	t_positions new_chain_end(unique_chains.n_elem);

	for (t_index i = 0; i < unique_chains.n_elem; ++i) {
		t_indices reads_in_chain(find(state.haplo == unique_chains(i)));
		new_haplo(reads_in_chain).fill(i);
		new_chain_start(i) = state.chain_start(unique_chains(i));
		new_chain_end(i) = state.chain_end(unique_chains(i));
	}

	state.haplo = new_haplo;
	state.chain_start = new_chain_start;
	state.chain_end = new_chain_end;
}



double haplotype::compute_delta_posterior(
		t_index const read_number,
		t_strand strand,
		t_haplochain const new_chain) const {

	TIMER_START
	DEBUG_ENTER

	t_haplochains new_haplo(state.haplo); //TODO this is alot of copying
	t_strands new_strands(state.read_strands); //TODO this is alot of copying
	t_positions new_starts(state.chain_start);
	t_positions new_ends(state.chain_end);
	double new_prior = state.prior;

	t_haplochain current_chain = state.haplo(read_number);

	if(current_chain == new_chain && strand == state.read_strands(read_number)) {
		return 0;
	}

	t_position region_start = data.reads_start_positions(read_number)-1;
	t_position region_end = data.reads_end_positions(read_number)+1;

	t_indices reads_in_region = find(data.reads_start_positions <= region_end && data.reads_end_positions >= region_start);


	//Create new haplo
	move_read(new_haplo, new_strands, new_starts, new_ends, new_prior, read_number, strand, new_chain);

	//Compute posterior

	return compute_delta(state.haplo, state.chain_start, state.chain_end, state.prior, state.read_strands,
			new_haplo, new_starts, new_ends, new_prior, new_strands,
			reads_in_region, region_start, region_end, current_chain, new_chain);
}

double haplotype::compute_delta_posterior(
		std::vector<t_indices> read_blocks,
		t_indices const& block,
		t_strand const strand,
		t_haplochain const new_chain) const {

	TIMER_START
	DEBUG_ENTER

	t_haplochains new_haplo(state.haplo);
	t_strands new_strands(state.read_strands); //TODO this is alot of copying
	t_positions new_starts(state.chain_start);
	t_positions new_ends(state.chain_end);
	double new_prior = state.prior;

	t_haplochain current_chain = state.haplo(block(0));

	if(current_chain == new_chain && strand == state.read_strands(block(0))) {
		return 0;
	}

	t_position region_start = min(data.reads_start_positions(block))-1; //FIXME
	t_position region_end = max(data.reads_end_positions(block))+1;

	t_indices reads_in_region = find(block_start_positions <= region_end && block_end_positions >= region_start);

	//create new haplo
	move_read_block(new_haplo, new_strands, new_starts, new_ends, new_prior, block, strand, new_chain);

	//Compute posterior
	return compute_delta(state.haplo, state.chain_start, state.chain_end, state.prior, state.read_strands,
			new_haplo, new_starts, new_ends, new_prior, new_strands,
			reads_in_region, region_start, region_end, current_chain, new_chain);
}


double haplotype::compute_delta(
		t_haplochains const& h,
		t_positions const& c_start,
		t_positions const& c_end,
		double prior,
		t_strands const& s,
		t_haplochains const& new_h,
		t_positions const& new_c_start,
		t_positions const& new_c_end,
		double new_prior,
		t_strands const& new_s,
		t_indices const& reads_in_region,
		t_position region_start,
		t_position region_end,
		t_haplochain current_chain,
		t_haplochain new_chain) const {

	TIMER_START
	DEBUG_ENTER

	double delta_posterior = 0;

	//prior
	delta_posterior -= prior_adjust*prior;
	delta_posterior += prior_adjust*new_prior;

	//current loglike
	delta_posterior -= compute_chain_loglike(h, s, current_chain, reads_in_region, region_start, region_end);

	//new loglike
	delta_posterior += compute_chain_loglike(new_h, new_s, current_chain, reads_in_region, region_start, region_end);

	if(new_chain != current_chain) {
		//current loglike
		delta_posterior -= compute_chain_loglike(h, s, new_chain, reads_in_region, region_start, region_end);

		//new loglike
		delta_posterior += compute_chain_loglike(new_h, new_s, new_chain, reads_in_region, region_start, region_end);

	}

	for(t_haplochain chain = max(h) + 1; chain <= max(new_h); ++chain) {

		if(chain == new_chain) {
			continue;
		}

		delta_posterior += compute_chain_loglike(new_h, new_s, chain, reads_in_region, region_start, region_end);
	}

#ifdef CHECK_POSTERIOR

	double delta = compute_posterior(new_h, new_c_start, new_c_end, new_s) - compute_posterior(h, c_start, c_end, s);
	if(fabs(delta_posterior - delta) > 1e-5) {

		cout << " ............ " << endl;
		cout << compute_chain_loglike(new_h, new_s, current_chain, reads_in_region, region_start, region_end)- compute_chain_loglike(h, s, current_chain, reads_in_region, region_start, region_end) << endl;
		cout << compute_chain_loglike(new_h, new_s, current_chain)- compute_chain_loglike(h, s, current_chain) << endl;

		throw std::runtime_error("compute_delta_posterior Error");
	}
#endif

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

	double logsum = 0;

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

	if(split_mode && is_DGCH(pos)) {

		// DGCH
		//rule out half methylated genotypes
		logsum(5, 0) = -1 * std::numeric_limits<float>::infinity();
		logsum(1, 4) = -1 * std::numeric_limits<float>::infinity();

	}

	else if(split_mode && is_HCGD(pos)) {

		logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
		logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	}

	else {

		//Non GpC or CpG context or normal mode

		//rule out half methylated genotypes
		logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
		logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	}


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

	t_indices reads_in_chain = find(state.haplo == chain);

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

		loglike_term.rows(read_start - chain_start, read_end - chain_start) += data.loglike_terms(*r, state.read_strands(*r));
	}

	for (t_position pos = chain_start; pos <= chain_end-1; ++pos) {
		compute_genotype(
				loglike_term.rows(pos - chain_start, pos - chain_start + 1),
				pos,
				genotype(pos-chain_start),
				genotype(pos+1-chain_start));
	}

	return genotype;

}

mat haplotype::compute_chain_loglike(t_haplochain const chain) const {

	t_indices reads_in_chain = find(state.haplo == chain);

	if (reads_in_chain.n_elem == 0) {
		throw std::runtime_error("compute_chain_loglike : chain empty");
	}

	//find start, end pos of chain
	t_positions starts = data.reads_start_positions(reads_in_chain);
	t_positions ends = data.reads_end_positions(reads_in_chain);

	t_position chain_start = min(starts);
	t_position chain_end = max(ends);

	mat loglike_term(chain_end - chain_start + 1, 6, fill::zeros);

	t_indices::const_iterator r = reads_in_chain.begin();
	for (; r != reads_in_chain.end(); ++r) {
		t_position const read_start = data.reads_start_positions(*r);
		t_position const read_end = data.reads_end_positions(*r);

		loglike_term.rows(read_start - chain_start, read_end - chain_start) += data.loglike_terms(*r, state.read_strands(*r));
	}

	return loglike_term;
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

	if(split_mode && is_DGCH(pos)) {

		// DGCH
		//rule out half methylated genotypes
		logsum(5, 0) = -1 * std::numeric_limits<float>::infinity();
		logsum(1, 4) = -1 * std::numeric_limits<float>::infinity();

	}

	else if(split_mode && is_HCGD(pos)) {

		logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
		logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	}

	else if(split_mode && pos > 0 && (is_DGCH(pos-1) || is_HCGD(pos-1))) {
		//Ignore
		return;
	}

	else {

		//Non GpC or CpG context or normal mode

		//rule out half methylated genotypes
		logsum(0, 5) = -1 * std::numeric_limits<float>::infinity();
		logsum(4, 1) = -1 * std::numeric_limits<float>::infinity();

	}

	t_index a;
	t_index b;
	logsum.max(a, b);

	g1 = static_cast<t_epi_base>(a);
	g2 = static_cast<t_epi_base>(b);

}


field<t_genotype> haplotype::genotypes() const {

	t_haplochains unique_chains = sort(unique(state.haplo));

	field<t_genotype> genotypes(unique_chains.n_elem);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			genotypes(i) = compute_chain_genotype(unique_chains(i));
		}

	return genotypes;
}

field<mat> haplotype::loglikes() const {

	t_haplochains unique_chains = sort(unique(state.haplo));

	field<mat> loglikes(unique_chains.n_elem);

		for (t_index i = 0; i < unique_chains.n_elem; ++i) {
			loglikes(i) = compute_chain_loglike(unique_chains(i));
		}

	return loglikes;
}

t_haplochains haplotype::compute_feasible_haplotypes(t_index id) const {

	DEBUG_ENTER
	TIMER_START

	t_haplochains feasible_haplotypes;

	//Compute for each chain the number of bases overlapping
	t_haplochains unique_chains = unique(state.haplo);

	for (t_index i = 0; i < unique_chains.n_elem; ++i) {

		if(use_read_blocks) {
			t_indices block = read_blocks[id];

			if(is_block_feasible(block, state.chain_start(unique_chains(i)), state.chain_end(unique_chains(i)))) {
				feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
				feasible_haplotypes(feasible_haplotypes.n_elem - 1) = unique_chains(i);
			}

		}

		else {

			t_position start = data.reads_start_positions(id);
			t_position end = data.reads_end_positions(id);

			if(is_read_feasible(start, end, state.chain_start(unique_chains(i)), state.chain_end(unique_chains(i)))) {
				feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
				feasible_haplotypes(feasible_haplotypes.n_elem - 1) = unique_chains(i);
			}

		}
	}

	//Add free haplo chain , i.e. a haplo chain with no reads
	feasible_haplotypes.resize(feasible_haplotypes.n_elem + 1);
	feasible_haplotypes(feasible_haplotypes.n_elem - 1) = max(state.haplo) + 1;

	return feasible_haplotypes;

}

double haplotype::compute_prior_chain(
		t_haplochains const& h,
		t_positions const& c_start,
		t_positions const& c_end,
		t_haplochain chain) const {

	TIMER_START
	DEBUG_ENTER

	uvec idx = find(h == chain);

	if(idx.n_elem == 0) {
		return 0;
	}

	double x = sum(data.reads_end_positions(idx) - data.reads_start_positions(idx) + 1);
	double l = c_end(chain) - c_start(chain)+1;

	return static_cast<double>(idx.n_elem) * (0.5*log(l)+2*log(x));

}

double haplotype::compute_posterior(
		t_haplochains const& h,
		t_positions const& c_start,
		t_positions const& c_end,
		t_strands const& strands) const {

	TIMER_START

	double posterior = 0;
	double p = 0;

	t_haplochains unique_chains = unique(h);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {

		//Compute prior
		p += compute_prior_chain(h, c_start, c_end, unique_chains(i));

		// Compute likelihood
		posterior += compute_chain_loglike(h, strands, unique_chains(i));

	}

	return posterior + prior_adjust*p;
}

double haplotype::compute_prior(
		t_haplochains const& h,
		t_positions const& c_start,
		t_positions const& c_end,
		t_strands const& strands) const {

	double p = 0;

	t_haplochains unique_chains = unique(h);
	for (t_index i = 0; i < unique_chains.n_elem; ++i) {

		//Compute prior
		p += compute_prior_chain(h, c_start, c_end, unique_chains(i));

	}

	return p;

}

double haplotype::posterior() const {
	return compute_posterior(state.haplo, state.chain_start, state.chain_end, state.read_strands);
}

#endif /* HAPLOTYPE_HPP_ */
