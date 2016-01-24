/*
 * alignment_data.hpp
 *
 *  Created on: Jan 26, 2014
 *      Author: martin
 */

#ifndef ALIGNMENT_DATA_HPP_
#define ALIGNMENT_DATA_HPP_

class alignment_data { //TODO better name for this class

public:

	t_count const n_reads;

	t_models const fwd_model;
	t_models const rev_model;

	t_models const fwd_DGCH_model;
	t_models const rev_DGCH_model;

	t_models const fwd_HCGD_model;
	t_models const rev_HCGD_model;

	t_models const fwd_C_G_model;
	t_models const rev_C_G_model;

	t_positions const reads_start_positions;
	t_positions const reads_end_positions;

	t_position const offset;
	t_position const sequence_length; //length of sequence

	t_counts const read_depth; //vector of size sequence length

	t_base_subsets const base_subsets; //vector of size (length of sequence)

	arma::field<arma::mat> const loglike_terms; //field of size (number of reads) x 2
	arma::field<arma::mat> const like_terms; //field of size (number of reads) x 2

	t_loglike_vector const log_pmax; //vector of size (number of reads)

	field<t_indices> const read_numbers; //Reads overlapping position
	t_indices const read_ids; //Reads ids
	std::vector<std::string> const read_names; //name of reads -- matching with read_id

	static bool is_HCGD(t_position pos, t_seq_bases const& ref) {
		return (pos > 0 && ref(pos-1) != 2 && ref(pos) == 1 && ref(pos+1) == 2 && ref(pos+2) != 1);
	}

	static bool is_DGCH(t_position pos, t_seq_bases const& ref) {
		return (pos > 0 && ref(pos-1) != 1 && ref(pos) == 2 && ref(pos+1) == 1 && ref(pos+2) != 2);
	}

    alignment_data(
    		AlgorithmConfiguration const& config,
    		std::vector<aligned_read> const& reads,
			std::string const refGenom_filename,
			std::string const refName,
			t_count hard_limit);

};

inline alignment_data::alignment_data(
		AlgorithmConfiguration const& config,
		std::vector<aligned_read> const& reads,
		std::string const refGenom_filename,
		std::string const refName,
		t_count hard_limit) :
        		n_reads(min(static_cast<t_count>(reads.size()), hard_limit)),
				fwd_model(config.fwd_model),
				rev_model(config.rev_model),
				fwd_DGCH_model(config.fwd_DGCH_model),
				rev_DGCH_model(config.rev_DGCH_model),
				fwd_HCGD_model(config.fwd_HCGD_model),
				rev_HCGD_model(config.rev_HCGD_model),
				fwd_C_G_model(config.fwd_C_G_model),
				rev_C_G_model(config.rev_C_G_model),
				offset(0),
				sequence_length(0),
				read_names(n_reads) {

	if(reads.empty()) {
		throw std::runtime_error("No reads");
	}

	t_positions starts(n_reads);
	t_positions ends(n_reads);
	std::vector<std::string> tmp_read_names(n_reads);
	t_indices tmp_read_ids(n_reads);

	t_count factor = max(static_cast<t_count>(reads.size())/n_reads,static_cast<t_count>(1));
	for (t_count i = 0; i < n_reads; ++i) {

		aligned_read read = reads[factor*i];

		starts(i) = read.position;
		ends(i) = read.position + read.length - 1;
		tmp_read_names[i] = read.name;
		tmp_read_ids(i) = factor*i;
	}

	//Set variables
	const_cast<t_position&>(this->offset) = starts.min();
	const_cast<t_position&>(this->sequence_length) = ends.max() - starts.min() + 1;

	const_cast<t_positions&>(this->reads_start_positions) = starts - offset;
	const_cast<t_positions&>(this->reads_end_positions) = ends - offset;

	//Load ref
    t_seq_bases ref = create_bases_vector(read_fasta(refGenom_filename, refName, offset, sequence_length+2));


	t_counts coverage(sequence_length, arma::fill::zeros);
	for (t_count i = 0; i < n_reads; ++i) {
		coverage.subvec(reads_start_positions(i), reads_end_positions(i)) += 1;
	}

	t_base_subsets subsets(sequence_length, arma::fill::zeros);
	arma::field<arma::mat> tmp_loglike_terms(n_reads, 2);
	arma::field<arma::mat> tmp_like_terms(n_reads, 2);

	arma::field<t_indices> tmp_read_numbers(sequence_length);

	t_loglike_vector tmp_pmax(n_reads);

	t_counts pos_read_count(sequence_length, arma::fill::zeros);

	for (t_count i = 0; i < n_reads; ++i) {

		aligned_read read = reads[factor*i];
		t_length L = read.length;

		//Check L <= max length
		if(L > fwd_model.n_elem) {
			throw std::runtime_error("Read exceeds max length");
		}

		tmp_loglike_terms(i, 0).set_size(L, 6);
		tmp_loglike_terms(i, 1).set_size(L, 6);

		tmp_like_terms(i, 0).set_size(L, 6);
		tmp_like_terms(i, 1).set_size(L, 6);

		for (t_position j = reads_start_positions(i);
				j <= reads_end_positions(i); ++j) {

			//read numbers
			if (tmp_read_numbers(j).is_empty()) {
				tmp_read_numbers(j).set_size(coverage(j));
			}

			tmp_read_numbers(j)(pos_read_count(j)) = i;
			++pos_read_count(j);

			//log sums
			t_position k = j - reads_start_positions(i);

			t_base base = read.bases(k);
			double epsilon = read.epsilon(k);

			if(base == 0) {
				//base = N
				tmp_loglike_terms(i, 0).row(k).fill(log(0.25));
				tmp_loglike_terms(i, 1).row(k).fill(log(0.25));

				tmp_like_terms(i, 0).row(k).fill(0.25);
				tmp_like_terms(i, 1).row(k).fill(0.25);

			} else {

				if(config.split_mode && j > 0 && (is_DGCH(j, ref) || (j > 1 && is_DGCH(j-1, ref)))) {
					//DGCH
					//fwd strand
					tmp_loglike_terms(i, 0).row(k) = log((1-4/static_cast<double>(3)*epsilon)*fwd_DGCH_model(k).row(base-1) + epsilon/static_cast<double>(3));

					//rev strand
					tmp_loglike_terms(i, 1).row(k) = log((1-4/static_cast<double>(3)*epsilon)*rev_DGCH_model(k).row(base-1) + epsilon/static_cast<double>(3));

					tmp_like_terms(i, 0).row(k) = (1-4/static_cast<double>(3)*epsilon)*fwd_DGCH_model(k).row(base-1) + epsilon/static_cast<double>(3);
					tmp_like_terms(i, 1).row(k) = (1-4/static_cast<double>(3)*epsilon)*rev_DGCH_model(k).row(base-1) + epsilon/static_cast<double>(3);

				}

				else if(config.split_mode && j > 0 && (is_HCGD(j, ref) || (j > 1 && is_HCGD(j-1, ref)))) {
					//HCGD
					tmp_loglike_terms(i, 0).row(k) = log((1-4/static_cast<double>(3)*epsilon)*fwd_HCGD_model(k).row(base-1) + epsilon/static_cast<double>(3));

					//rev strand
					tmp_loglike_terms(i, 1).row(k) = log((1-4/static_cast<double>(3)*epsilon)*rev_HCGD_model(k).row(base-1) + epsilon/static_cast<double>(3));

					tmp_like_terms(i, 0).row(k) = (1-4/static_cast<double>(3)*epsilon)*fwd_HCGD_model(k).row(base-1) + epsilon/static_cast<double>(3);
					tmp_like_terms(i, 1).row(k) = (1-4/static_cast<double>(3)*epsilon)*rev_HCGD_model(k).row(base-1) + epsilon/static_cast<double>(3);
				}

				else if(config.split_mode && (ref(j) == 1 || ref(j) == 2)) {
					//GC CG or single C G
					tmp_loglike_terms(i, 0).row(k) = log((1-4/static_cast<double>(3)*epsilon)*fwd_C_G_model(k).row(base-1) + epsilon/static_cast<double>(3));

					//rev strand
					tmp_loglike_terms(i, 1).row(k) = log((1-4/static_cast<double>(3)*epsilon)*rev_C_G_model(k).row(base-1) + epsilon/static_cast<double>(3));

					tmp_like_terms(i, 0).row(k) = (1-4/static_cast<double>(3)*epsilon)*fwd_C_G_model(k).row(base-1) + epsilon/static_cast<double>(3);
					tmp_like_terms(i, 1).row(k) = (1-4/static_cast<double>(3)*epsilon)*rev_C_G_model(k).row(base-1) + epsilon/static_cast<double>(3);
				}

				else {
					//Non GpC context or normal mode

					//fwd strand
					tmp_loglike_terms(i, 0).row(k) = log((1-4/static_cast<double>(3)*epsilon)*fwd_model(k).row(base-1) + epsilon/static_cast<double>(3));

					//rev strand
					tmp_loglike_terms(i, 1).row(k) = log((1-4/static_cast<double>(3)*epsilon)*rev_model(k).row(base-1) + epsilon/static_cast<double>(3));

					tmp_like_terms(i, 0).row(k) = (1-4/static_cast<double>(3)*epsilon)*fwd_model(k).row(base-1) + epsilon/static_cast<double>(3);
					tmp_like_terms(i, 1).row(k) = (1-4/static_cast<double>(3)*epsilon)*rev_model(k).row(base-1) + epsilon/static_cast<double>(3);

				}

			}

			//Position marks
			subsets(j) = base_union(subsets(j), base);
		}

		tmp_pmax(i) = max(as_scalar(mean(max(tmp_loglike_terms(i,0), 1))), as_scalar(mean(max(tmp_loglike_terms(i,1), 1))));
	}

	//Set variables
	const_cast<t_base_subsets&>(this->base_subsets) = subsets;
	const_cast<arma::field<arma::mat>&>(this->loglike_terms) = tmp_loglike_terms;
	const_cast<arma::field<arma::mat>&>(this->like_terms) = tmp_like_terms;
	const_cast<field<t_indices>&>(this->read_numbers) = tmp_read_numbers;
	const_cast<t_indices&>(this->read_ids) = tmp_read_ids;
	const_cast<std::vector<std::string>&>(this->read_names) = tmp_read_names;
	const_cast<t_loglike_vector&>(this->log_pmax) = tmp_pmax;
	const_cast<t_counts&>(this->read_depth) = coverage;
}




#endif /* ALIGNMENT_DATA_HPP_ */
