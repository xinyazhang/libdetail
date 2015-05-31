#ifndef KDTREE_H
#define KDTREE_H

#include <algorithm>
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include "bbox.h"
#include "ray.h"
#include "../pref.h"

using std::cerr;
using std::endl;

class Geometry;

const int KDTREE_SPLIT_THRESH = 10;
const double KDTREE_IMBA_THRESH = 0.8;
const int KDTREE_DIVIDER_SEARCH_THRESH = 5;
const double KDTREE_DIVIDER_BALANCE_THRESH = 0.1;
const double KDTREE_OVERKILL_THRESH = 8.0;
//const int KDTREE_DEPTH = 12;

extern bool debugMode;
extern const void* debugHitObj;
void debugDummy();

struct isect;
static constexpr const char axis_name[] = "XYZ";

struct line_endpoint {
	double coord;
	int dlobj;
	int drobj;

#if 0
	line_endpoint(double x, int y, int z)
		:coord(x), dlobj(y), drobj(z)
	{
	}
#endif
};

inline bool operator < (const line_endpoint& lhs, const line_endpoint& rhs)
{
	return lhs.coord < rhs.coord;
}

/*
 * Concept of G:
 * 	const BoundingBox& getBoundingBox() const;
 *
 */
template<typename G>
class KdTreeNode {
private:
	double axis(const Vec3d& vec)
	{
		return vec[axis_id_];
	}
	double bbmin() { return axis(bounds().getMin()); }
	double bbmax() { return axis(bounds().getMax()); }
public:
	KdTreeNode(int depth = 0, int initial_axis = 0)
		:depth_(depth), axis_id_(initial_axis)
	{
		divider_ = NAN;
	}
	size_t size() const { return objs_.size(); }
	void add(G* g)
	{
		objs_.emplace_back(g);
		bounds_.merge(g->getBoundingBox());
	}

	void split(int budget)
	{
		if (size() <= KDTREE_SPLIT_THRESH) {
			return ;
		}

		if (budget > 0 && search_divider()) {
			children_[0].reset(new KdTreeNode<G>(depth_ + 1, mawaru_axis()));
			children_[1].reset(new KdTreeNode<G>(depth_ + 1, mawaru_axis()));
			double actualoverkill = add_objs_to_children();
			if (actualoverkill <= KDTREE_OVERKILL_THRESH) {
				for(int i = 0; i < 2; i++) {
					children_[i]->split(budget - 1);
				}
			} else
				divider_ = NAN;
		} else
			divider_ = NAN;
	}
	
	double add_objs_to_children()
	{
		size_t added = 0;
		for(const auto& g : objs_) {
			const BoundingBox& bb = g->getBoundingBox();
			bool oops = true;
			if (axis(bb.getMin()) <= divider_) { // "Left" 
				children_[0]->add(g);
				oops = false;
				added++;
			}
			if (axis(bb.getMax()) >= divider_) { // "Right"
				children_[1]->add(g);
				oops = false;
				added++;
			}
#ifndef NDEBUG
			if (oops)
				for (int i = 0 ; i < 120; i++) cerr << "##############################OOPS" << endl;
#endif
		}
		int nlobj = children_[0]->size();
		int nrobj = children_[1]->size();
#ifndef NDEBUG
		for(int i = 0; i < depth_; i++) cerr << "  ";
		cerr << objs_.size() << " was ACTUALLY divided into "
			<< nlobj<<", "<< nrobj
			<<" on Axis " << axis_name[axis_id_];
#endif
		double imba = check_imba(divider_, nlobj, nrobj);
#ifndef NDEBUG
		cerr << " with IMBA: " << imba;
#endif
		if (imba >= KDTREE_OVERKILL_THRESH) {
#ifndef NDEBUG
			cerr << "\t IMBA, stopped here" << endl;
#endif
			return KDTREE_OVERKILL_THRESH + 2.0;
		}
#ifndef NDEBUG
		cerr << endl;
#endif
		return double(added)/objs_.size();
	}

	const BoundingBox& bounds() const { return bounds_; }

	bool simple_search(ray& r, isect& i)
	{
		bool have_one = false;
		for(auto g : objs_) {
			isect cur;
			if (g->intersect(r, cur)) {
				if (!have_one || (cur.t < i.t)) {
					i = cur;
					have_one = true;
				}
			}
		}
		if(!have_one)
			i.setT(1e308);
		return have_one;
	}


	bool intersect(ray& r, isect& is)
	{
		if (size() == 0)
			return false;
		if (std::isnan(divider_))
			return simple_search(r, is);

		bool hitbox[2];
		bool both = false;
		double tmin[2], tmax[2];
		int first = 0, second = 1; // child indices to check
		for(int i = 0; i < 2; i++) {
			if (!children_[i]) {
				hitbox[i] = false;
				continue;
			}
			const BoundingBox& bb = children_[i]->bounds();
			hitbox[i] = bb.intersect(r, tmin[i], tmax[i]);

			if (tmax[i] > RAY_EPSILON && tmin[i] < RAY_EPSILON)
				tmin[i] = tmax[i];
		}
		bool ret = false; 
		for (int child = 0; child < 2; child++) {
			isect probe;
			if (ret && is.t < tmin[child])
				continue;
			if (intersect_child(child, r, hitbox[child], probe)) {
				if (!ret || (probe.t < is.t)) {
					is = probe;
					ret = true;
				}
			}
		}
		return ret;

		return ret;
	}
	bool hunter_obj(const void* target)
	{
		if (!size())
			return false;
		bool ret = has(target);
		if (!ret)
			return false;
		report_self(); cerr << endl;

		for(int i = 0; i < 2; i++) {
			if (!children_[i])
				continue;
			if (!children_[i]->hunter_obj(target))
				continue;
		}
		return ret;
	}
	bool has(const void* target)
	{
		bool ret = false;
		for(const auto& ptr : objs_)
			if (ptr == target)
				ret = true;
		return ret;
	}

	void report_self()
	{
#ifndef NDEBUG
		const BoundingBox& bb = bounds();
		for(int i = 0; i < depth_; i++) cerr << "  ";
		cerr << "Size: " << size() << ", Axis " << axis_name[axis_id_];
		cerr << " divided by " << divider_ << "\t";
		cerr << "Find object " << debugHitObj
			<< " with BB " << bb.getMin()
			<< " --- "
			<< bb.getMax();
#endif
	}
private:
	std::vector<G*> objs_;
	std::unique_ptr<KdTreeNode<G>> children_[2];
	BoundingBox bounds_; // only includes objecs
	int axis_id_;
	double divider_;
	int depth_;

	bool intersect_child(int idx, ray& r, bool hitbox, isect& probe)
	{
		static const char axisign[] = "-+";
		if (!hitbox)
			return false;
		bool ret;
		ret = children_[idx]->intersect(r, probe);

		return ret;
	}

	void calc_nobjs(double middle,
			int& nlobj,
			int& nrobj)
	{
		nlobj = 0;
		nrobj = 0;
		for(const auto& g : objs_) {
			const BoundingBox& bb = g->getBoundingBox();
			if (axis(bb.getMin()) <= middle) { // "Left" 
				nlobj++;
			}
			if (axis(bb.getMax()) >= middle) { // "Right"
				nrobj++;
			}
		}
	}

	bool adjust_divider(double middle,
			int& nlobj, int& nrobj)
	{
		double lbound = axis(bounds().getMin());
		double rbound = axis(bounds().getMax());
		ssize_t imbasize = (ssize_t)objs_.size() * KDTREE_DIVIDER_BALANCE_THRESH;
		double best_imba = KDTREE_OVERKILL_THRESH * 2;
		double best_middle = middle;
		double imba;
		int niter = 0;
		do {
			calc_nobjs(middle, nlobj, nrobj);
			imba = check_imba(middle, nlobj, nrobj);
			if (imba < best_imba) {
				best_middle = middle;
				best_imba = imba;
			}
#ifndef NDEBUG
#if 0
			cerr << "\t\t\tAdjust divider on Axis "<< axis_name[axis_id_]
				<< ", middle: " << middle << endl;
#endif
#endif
			if (nlobj - nrobj > imbasize) {
				rbound = middle;
			} else if (nrobj - nlobj > imbasize) {
				lbound = middle;
			} else
				break; // balanced
			middle = (lbound + rbound) / 2;
			niter++;
		} while (niter < KDTREE_DIVIDER_SEARCH_THRESH);
		imba = check_imba(middle, nlobj, nrobj);
		if (imba < best_imba) {
			best_middle = middle;
			best_imba = imba;
		}
		if (best_imba > KDTREE_OVERKILL_THRESH)
			return false;
		divider_ = best_middle;
		return true;
	}

	bool set_divider(double& overkill)
	{
		std::vector<line_endpoint> axes(objs_.size() * 2);
		double middle = 0.5*(bbmin() + bbmax());

		
		bool ret = true;
		overkill = 1.0;
		int nlobj, nrobj;
		ret = adjust_divider(middle, nlobj, nrobj);
		overkill = check_imba(middle, nlobj, nrobj);
		return ret;
	}

	double check_imba(double middle, int nlobj, int nrobj)
	{
		double span = bbmax() - bbmin();
		double rate = (middle - bbmin())/span;
		int larger = std::max(nlobj, nrobj);
		if (larger == size())
			return KDTREE_OVERKILL_THRESH * 2;
		double x = size()/double(size()-larger);
		return std::log(x) * (double(nlobj + nrobj)/size());
	}

	bool search_divider()
	{
		int best_axis = axis_id_;
		double bestoverkill = KDTREE_OVERKILL_THRESH + 2.0;
		bool ret = false;
		double overkill;

		for(int i = 0; i < 3; i++) {
			if (set_divider(overkill)) {
				ret = true;
				if (bestoverkill > overkill) {
					bestoverkill = overkill;
					best_axis = axis_id_;
				}
			}
			axis_id_ = mawaru_axis();
		}
		axis_id_ = best_axis;
		return set_divider(overkill);
	}

	int mawaru_axis() const
	{
		return (axis_id_ + 1) % 3;
	}
};

/*
 * Handle special cases before entering the Tree.
 * 	e.g. Check whether the ray intersects with the outmost bounding box
 * 	     , which isn't done by KdTreeNode since the node only checks its
 * 	     children's bounding box.
 */
template<typename G> // G: Geometry
class KdTree {
public:
	void add(G* g)
	{
		root_.add(g);
		objects.emplace_back(g);
	}

	void split()
	{
		root_.split(pref::getKDTreeDepth(root_.size()));
	}

	bool intersect(ray& r, isect& i)
	{
		bool have_one = false;
		double tmin = 0.0;
		double tmax = 0.0;
		if (debugMode || !pref::kdtree_enabled) {
			have_one = false;
			debugHitObj = nullptr;
			for(auto j : objects) {
				isect cur;
				if( j->intersect(r, cur) ) {
					if(!have_one || (cur.t < i.t)) {
						i = cur;
						have_one = true;
					}
				}
			}
			if(!have_one) {
				i.setT(1000.0);
			} else if (debugMode) {
#ifndef NDEBUG
				const BoundingBox& bb = i.obj->getBoundingBox();
				cerr << "Reference object's bounding box "
					<< bb.getMin()
					<< " --- "
					<< bb.getMax()
					<< endl;
#endif
			}
		}
		if (pref::kdtree_enabled) {
			const void* target = i.obj;
			debugHitObj = i.obj;
			if (root_.bounds().intersect(r, tmin, tmax))
				have_one = root_.intersect(r, i);
			else
				have_one = false;
#ifndef NDEBUG
			if (debugMode && i.obj && !have_one) {
				cerr << "### Answer Mismatch for Ray (" << r.getPosition()
					<<") => ("<< r.getDirection() << ")" << endl;
			}
#endif
		}
		return have_one;
	}
	bool hunter_obj(const void* target)
	{
#ifndef NDEBUG
		cerr << "=============== Object Hunter ===============" << endl;
		root_.hunter_obj(target);
#endif
	}
private:
	KdTreeNode<G> root_;
	std::vector<G*> objects;

	void dummy() {}
};

#endif
