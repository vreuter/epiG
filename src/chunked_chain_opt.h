/*
 * openmp_optimizer.hpp
 *
 *  Created on: Nov 15, 2013
 *      Author: martin
 */

#ifndef CHUNK_OPTIMIZER_HPP_
#define CHUNK_OPTIMIZER_HPP_

using namespace arma;

class chunk_haplo_chain_optimizer {

public:
	t_count const number_of_chunks;

	field<std::string> const refNames;

	t_position const offset;

	t_count const max_threads;

	field<AlgorithmConfiguration> const configs;

private:

	std::string bam_file;

	field<field<t_indices> > read_numbers; // one for each chunk
	field<t_indices> read_ids;
	field<std::vector<std::string> > read_names;

	field<t_haplochains> haplotypes;
	field<t_strands> strands;
	field<field<Col<t_epi_base> > > genotypes;
	field<field<mat> > loglikes;
	field<t_positions> chain_starts;
	field<t_positions> chain_ends;

	t_positions chunk_start_pos;
	t_positions chunk_end_pos;

public:

	chunk_haplo_chain_optimizer(
			std::string const& bam_file,
			field<std::string> const& refNames,
			t_positions const& chunk_start_positions,
			t_positions const& chunk_end_positions,
			t_count max_threads,
			field<AlgorithmConfiguration> const& configs);

	void run();

	field<t_haplochains> get_haplotypes() const {
		return haplotypes;
	}

	field<t_strands> get_strands() const {
		return strands;
	}

	field<field<Col<t_epi_base> > > get_genotypes() const {
		return genotypes;
	}

	field<field<mat> > get_loglikes() const {
		return loglikes;
	}

	field<t_positions> get_chain_start() const {
		return chain_starts;
	}

	field<t_positions> get_chain_end() const {
		return chain_ends;
	}

	field< field<t_indices> > get_readNumbers() const {
		return read_numbers;
	}

	field<t_indices> get_readIDs() const {
		return read_ids;
	}

	field<std::vector<std::string> > get_readNames() const {
		return read_names;
	}


	t_positions get_chunks_start() const {
		return chunk_start_pos;
	}

	t_positions get_chunks_end() const {
		return chunk_end_pos;
	}
};

chunk_haplo_chain_optimizer::chunk_haplo_chain_optimizer(
		std::string const& bam_file,
		field<std::string> const& refNames,
		t_positions const& chunk_start_positions,
		t_positions const& chunk_end_positions,
		t_count max_threads,
		field<AlgorithmConfiguration> const& configs) :
				number_of_chunks(chunk_start_positions.n_elem),
				refNames(refNames),
				offset(chunk_start_positions.min()),
				max_threads(max_threads),
				configs(configs),
				bam_file(bam_file),
				read_numbers(number_of_chunks),
				read_ids(number_of_chunks),
				read_names(number_of_chunks),
				haplotypes(number_of_chunks),
				strands(number_of_chunks),
				genotypes(number_of_chunks),
				loglikes(number_of_chunks),
				chain_starts(number_of_chunks),
				chain_ends(number_of_chunks),
				chunk_start_pos(chunk_start_positions),
				chunk_end_pos(chunk_end_positions) {

	//TODO domain check consistency of chunk positions
	//refnames length =n_chunks
	// length configs = n_chunks
}

void chunk_haplo_chain_optimizer::run() {

	// create progress monitor
	Progress p(number_of_chunks, configs(0).verbose);

	//Warnings
	omp_rwarn warnings;

#ifdef EPIG_USE_OPENMP
	omp_set_num_threads(max_threads);
#pragma omp parallel for schedule(dynamic)
#endif
	for (t_count i = 0; i < number_of_chunks; ++i) {

		const AlgorithmConfiguration config(configs(i));
		long unsigned int n_reads_hard_limit = config.reads_hard_limit;

		if (!p.is_aborted()) {

			bamReader reader(bam_file, refNames(i), chunk_start_pos(i),
					chunk_end_pos(i));

#ifdef EPIG_USE_OPENMP
#pragma omp critical
#endif
			{
				//TODO does this need to be in a critical region ??
				//TODO use limit when fetching reads

				reader.fetch(); //Load reads

				//remove low quality reads
				reader.remove_low_quality(config.quality_threshold);

				//Check that chunk is non empty
				if (reader.get_reads().empty()) {
					read_ids(i).set_size(0);
					haplotypes(i).set_size(0);
				}

			}

			if (reader.get_reads().empty()) {
				p.increment();
				continue;
			}

			if (reader.get_reads().size() > n_reads_hard_limit) {
				std::ostringstream msg;
				msg << reader.get_reads().size()
						<< " reads are covering chunk (position: "
						<< chunk_start_pos(i) << " : " << chunk_end_pos(i)
						<< ") this exceeds the hard limit. Only "
						<< n_reads_hard_limit << " reads are used.";
				warnings.add(msg.str());
			}
	        //Load and prepare data
	        std::vector<aligned_read> const& reads = reader.get_reads();
	        alignment_data data(config, reads, refNames(i), n_reads_hard_limit);

			//Load ref
			t_seq_bases ref;
			if(config.use_ref) {
				ref = create_bases_vector(read_fasta(config.ref_filename, config.ref_offset, refNames(i), data.offset, data.sequence_length+2));

	      //TODO append 0's (unknowns) + error -> warning
	      if (ref.n_elem != static_cast<unsigned int>(data.sequence_length+2)) {
					throw std::runtime_error("Problem with refGenom"); //TODO error msg
				}
			}	

			else {
				ref.resize(data.sequence_length+2);
				ref.zeros();
			}

	    //Load alt
			t_seq_bases alt;
			if(config.use_alt) {

				alt = create_bases_vector(read_fasta(config.alt_filename, config.alt_offset,  refNames(i), data.offset, data.sequence_length+2));

	      //TODO append 0's (no snp) + error -> warning
	      if (alt.n_elem != static_cast<unsigned int>(data.sequence_length+2)) {
	      	std::ostringstream msg;
	      	msg << "Problem with altGenom"; //TODO warning msg
	      	warnings.add(msg.str());

        	int old_n_elem = alt.n_elem;
        	alt.resize(data.sequence_length+2);
        	alt.subvec(old_n_elem, alt.n_elem-1).zeros();
				}
			}

			else {
				alt.resize(data.sequence_length+2);
				alt.zeros();
			}

			//Init haplo optimizer
			haplo_chain_optimizer opt = haplo_chain_optimizer::create(config, data, ref, alt, warnings, config.use_paired_reads);

			opt.run(p);

#ifdef EPIG_USE_OPENMP
#pragma omp critical
#endif
			{

				chunk_start_pos(i) = max(chunk_start_pos(i),
						data.reads_start_positions.min() + data.offset);
				chunk_end_pos(i) = min(chunk_end_pos(i),
						data.reads_end_positions.max() + data.offset);

				read_numbers(i) = data.read_numbers.subfield(
						chunk_start_pos(i) - data.offset, 0,
						chunk_end_pos(i) - data.offset, 0);

				read_ids(i) = data.read_ids;

				read_names(i) = data.read_names;

				haplotypes(i) = opt.get_haplochains();

				strands(i) = opt.get_read_strands();

				genotypes(i) = opt.get_chain_genotypes();

				loglikes(i) = opt.get_chain_loglikes();

				chain_starts(i) = opt.haplo_chain_start() + data.offset ;

				chain_ends(i) = opt.haplo_chain_end() + data.offset;

			}

			//TODO remove
			//cout << i << " done" << endl;

			p.increment();
		}
	}

}


chunk_haplo_chain_optimizer create_base_chunk_optimizer(
		std::string const& bam_file,
		std::string const& refName,
		t_position start_position,
		t_position end_position,
		t_count max_threads, t_position chunk_size,
		field<AlgorithmConfiguration> const& configs) {

	t_count number_of_chunks = ceil(
			(end_position - start_position + 1)
					/ static_cast<double>(chunk_size));

	t_positions chunk_start_pos(number_of_chunks);
	t_positions chunk_end_pos(number_of_chunks);

	//Compute chunk positions
	chunk_start_pos(0) = start_position;
	chunk_end_pos(0) = chunk_start_pos(0) + chunk_size - 1;

	for (t_count i = 1; i < number_of_chunks; ++i) {
		chunk_start_pos(i) = chunk_end_pos(i - 1) + 1;
		chunk_end_pos(i) = chunk_start_pos(i) + chunk_size - 1;
	}

	//Correct last end position
	chunk_end_pos(number_of_chunks - 1) = end_position;

	//Ref names
	field<std::string> refNames(number_of_chunks);
	for (t_count i = 0; i < number_of_chunks; ++i) {
		refNames(i) = refName;
	}

	return chunk_haplo_chain_optimizer(bam_file, refNames, chunk_start_pos, chunk_end_pos, max_threads, configs);
}


chunk_haplo_chain_optimizer create_base_chunk_optimizer(
		std::string const& bam_file,
		std::string const& refName,
		t_position start_position,
		t_position end_position,
		t_count max_threads,
		t_position chunk_size,
		AlgorithmConfiguration const & config) {

	field<AlgorithmConfiguration> configs(1);
	configs(0) = config;

	return create_base_chunk_optimizer(
			bam_file,
			refName,
			start_position,
			end_position,
			max_threads,
			chunk_size,
			configs);
}

#endif /* CHUNK_OPTIMIZER_HPP_ */
