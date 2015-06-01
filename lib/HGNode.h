#ifndef HGNODE_H
#define HGNODE_H

#include <memory>
#include <vector>

namespace hgeo {

using namespace std;

/*
 * HGNode: container of geometries
 *
 *	Space: the mathematical space that the geometry belongs to
 *	Geometry: type of the objects in this space
 *	Partitioner: class that split the Space into subspaces.
 *
 * User can inherit this class to implement algorithms over those geometries
 */
template<typename Space, typename Geometry, typename Partitioner>
class HGNode {
	typedef HGNode<Space, Geometry, Partitioner> self_type;
public:
	HGNode(const Space& = Space());
	~HGNode();

	// Geometry manipulation
	void add(Geometry* geo) { geos_.emplace_back(geo); }
	const vector<const Geometry*>& get_geos() const { return geos_; }
	vector<Geometry*> get_geos() { return geos_; }

	// Child manipulation
	shared_ptr<self_type> child_at(size_t idx);
	bool has_child() const { return !children_.empty(); }
	void redistribute();

protected:
	bool create_child_at(size_t idx, Partitioner* ppart = nullptr);

	Space space_;
	vector<Geometry*> geos_;
	vector<shared_ptr<self_type>> children_;
};

#define HGNODE_CLASS HGNode<Space,Geometry,Partitioner>
#define HGNODE_SCOPE(ret_type) template<typename Space, typename Geometry, typename Partitioner> \
ret_type HGNODE_CLASS::

HGNODE_SCOPE(void)
create_child_at(size_t idx, Partitioner* ppart = nullptr)
{
	if (children_.size() > idx && children_[idx])
		return false; // Do not recreate child
	if (!ppart) {
		Partitioner part(this->space_, geos_);
		create_child_at(idx, &part);
	}
	if (!ppart->set_subspace(idx))
		return false;
	auto child = shared_ptr<self_type>(new HGNode(ppart->get_space()));
	Geometry* geo;
	while (nullptr != (geo = part.get_next())) {
		child->add(geo);
	}
	if (child->size() == 0)
		return false;
	if (children_.size() <= idx)
		children_.resize(idx+1);
	children_[idx] = child;
	return true;
}

HGNODE_SCOPE(void)
redistribute()
{
	children_.clear();
	Partitioner part(this->space_, geos_);
	for (size_t i = 0; i < part.get_nsubspace(); i++)
		create_child_at(idx, &part);
}

HGNODE_SCOPE(shared_ptr<HGNODE_CLASS>)
child_at(size_t idx)
{
	create_child_at(idx);
	if (children_.size() > idx)
		return children_[idx];
	return shared_ptr<self_type>(); // null
}

#undef HGNODE_SCOPE
#undef HGNODE_CLASS

};

#endif
