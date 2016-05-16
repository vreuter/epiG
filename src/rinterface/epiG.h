#ifndef EPIG_HAPLO_R_INTERFACE_HPP_
#define EPIG_HAPLO_R_INTERFACE_HPP_

extern "C" {

SEXP r_epiG_haplo_fit_filename(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_max_threads, SEXP r_max_chunk_size, SEXP r_config);
SEXP r_epiG_haplo_fit_filename_chunks(SEXP r_filename, SEXP r_refName, SEXP r_chunks_start, SEXP r_chunks_end, SEXP r_max_threads, SEXP r_config);

SEXP r_epiG_compute_chunk_positions(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_chunk_size);

}


SEXP epiG_compute_chunk_positions(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_chunk_size) {

    const std::string filename = get_value<std::string>(r_filename);
    const std::string refName = get_value<std::string>(r_refName);
    const t_position start = get_value<t_position>(r_start);
    const t_position end = get_value<t_position>(r_end);
    const unsigned int chunk_size = get_value<unsigned int>(r_chunk_size);

    bamReader reader(filename, refName, start, end);
    reader.fetch_reads_info();
    std::vector<read_info> const& infos = reader.get_infos();

    //t_positions start_pos(infos.size());
    t_positions end_pos(infos.size());

    unsigned int n_reads = infos.size();

    for(unsigned int i = 0; i < n_reads; ++i) {

        read_info info = infos[i];

      //  start_pos(i) = info.position;
        end_pos(i) = info.position + info.length - 1;
    }

    std::vector<t_position> pos;
    pos.push_back(0);

    unsigned int end_sum = 0;

    for(unsigned int i = 0; i < n_reads; ++i) {

        if((i+1 - end_sum) >= chunk_size) {
                pos.push_back(end_pos(i)+1);
                end_sum = sum(end_pos.subvec(0,i-1) < pos.back());
        }
    }

    //Correct ends
    pos[0] = start;
    pos.back() = end;

    return rObject(pos);
}

SEXP r_epiG_compute_chunk_positions(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_chunk_size) {

    try {

        return epiG_compute_chunk_positions(r_filename, r_refName, r_start, r_end, r_chunk_size);

        //Catch unhandled exceptions

    } catch (std::exception & e) {

        if (e.what() != NULL) {
            report_error(e.what());
        }

        else {
            report_error("Unknown error");
        }

    } catch (...) {
        report_error("Unknown error");
    }

    return R_NilValue; //Avoid compiler warnings
}


SEXP epiG_haplo_fit_filename(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_max_threads, SEXP r_max_chunk_size, SEXP r_config) {

	const std::string filename = get_value<std::string>(r_filename);
	const std::string refName = get_value<std::string>(r_refName);
	const t_position start = get_value<t_position>(r_start);
	const t_position end = get_value<t_position>(r_end);
	const t_count max_threads = get_value<t_count>(r_max_threads);
	const t_count max_chunk_size = get_value<t_count>(r_max_chunk_size);

	const AlgorithmConfiguration config = get_value<AlgorithmConfiguration>(r_config);

	chunk_haplo_chain_optimizer opt = create_base_chunk_optimizer(filename, refName, start, end, max_threads, max_chunk_size, config);

	opt.run();

    rList res;

	res.attach(rObject(opt.offset), "offset");
	res.attach(rObject(opt.number_of_chunks), "number_of_chunks");
	res.attach(rObject(opt.get_chunks_start()), "chunks_start");
	res.attach(rObject(opt.get_chunks_end()), "chunks_end");
	res.attach(rObject(opt.get_readNumbers()), "read_ids");
	res.attach(rObject(opt.get_readIDs()), "read_unique_chunk_ids");
	res.attach(rObject(opt.get_readNames()), "read_names");
	res.attach(rObject(opt.get_haplotypes()), "haplotypes");
	res.attach(rObject(opt.get_strands()), "strands");
	res.attach(rObject(opt.get_genotypes()), "genotypes");
	res.attach(rObject(opt.get_loglikes()), "loglikes");
	res.attach(rObject(opt.get_chain_start()), "chain_start");
	res.attach(rObject(opt.get_chain_end()), "chain_end");

    return rObject(res);
}

SEXP r_epiG_haplo_fit_filename(SEXP r_filename, SEXP r_refName, SEXP r_start, SEXP r_end, SEXP r_max_threads, SEXP r_max_chunk_size, SEXP r_config) {

	try {

		return epiG_haplo_fit_filename(r_filename, r_refName, r_start, r_end, r_max_threads, r_max_chunk_size, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}

SEXP epiG_haplo_fit_filename_chunks(SEXP r_filename, SEXP r_refNames, SEXP r_chunks_start, SEXP r_chunks_end, SEXP r_max_threads, SEXP r_configs) {

	const std::string filename = get_value<std::string>(r_filename);
	const field<std::string> refNames = get_field<std::string>(r_refNames);
	const t_positions chunks_start = get_value<t_positions>(r_chunks_start);
	const t_positions chunks_end = get_value<t_positions>(r_chunks_end);
	const t_count max_threads = get_value<t_count>(r_max_threads);

	const field<AlgorithmConfiguration>  configs = get_field<AlgorithmConfiguration>(r_configs);

	chunk_haplo_chain_optimizer opt(filename, refNames, chunks_start, chunks_end, max_threads, configs);

	opt.run();

    rList res;

	res.attach(rObject(opt.offset), "offset");
	res.attach(rObject(opt.number_of_chunks), "number_of_chunks");
	res.attach(rObject(opt.get_chunks_start()), "chunks_start");
	res.attach(rObject(opt.get_chunks_end()), "chunks_end");
	res.attach(rObject(opt.get_readNumbers()), "read_ids");
	res.attach(rObject(opt.get_readIDs()), "read_unique_chunk_ids");
	res.attach(rObject(opt.get_haplotypes()), "haplotypes");
	res.attach(rObject(opt.get_strands()), "strands");
	res.attach(rObject(opt.get_genotypes()), "genotypes");
	res.attach(rObject(opt.get_loglikes()), "loglikes");
	res.attach(rObject(opt.get_chain_start()), "chain_start");
	res.attach(rObject(opt.get_chain_end()), "chain_end");

    return rObject(res);
}

SEXP r_epiG_haplo_fit_filename_chunks(SEXP r_filename, SEXP r_refNames, SEXP r_chunks_start, SEXP r_chunks_end, SEXP r_max_threads, SEXP r_config) {

	try {

		return epiG_haplo_fit_filename_chunks(r_filename, r_refNames, r_chunks_start, r_chunks_end, r_max_threads, r_config);

		//Catch unhandled exceptions

	} catch (std::exception & e) {

		if (e.what() != NULL) {
			report_error(e.what());
		}

		else {
			report_error("Unknown error");
		}

	} catch (...) {
		report_error("Unknown error");
	}

	return R_NilValue; //Avoid compiler warnings
}


#endif /* EPIG_HAPLO_R_INTERFACE_HPP_ */
