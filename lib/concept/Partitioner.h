#ifndef CONCEPT_PARTITIONER_H
#define CONCEPT_PARTITIONER_H
/*
 * This file describes the concept of class Partitioner
 * Users should define their own classes which provide the same interface.
 * DO NOT INHERIT this class directly.
 */

template<typename Space, typename Geometry>
class Partitioner {
public:
	Partitioner(const Space& parent_space,
			std::vector<Geometry*> geos_in_parent_space);
	size_t get_nsubspace() const;
	void set_subspace(size_t idx);

	Space get_space();
	Geometry* get_next();
};

#endif
