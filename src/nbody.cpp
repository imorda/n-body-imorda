#include <cassert>
#include <cmath>
#include <fstream>
#include <nbody.h>

Body::Body(const Cartesian & m_pos, const Cartesian & m_speed, const double m_mass, const std::string & mName)
    : m_pos(m_pos)
    , m_speed(m_speed)
    , m_mass(m_mass)
    , m_name(mName)
{
}
Body::Body()
    : m_pos{0, 0}
    , m_speed{0, 0}
    , m_mass(0)
{
}
double Body::distance(const Body & b) const
{
    return hypot(m_pos.x - b.m_pos.x, m_pos.y - b.m_pos.y);
}
void Body::add_force(const Body & b)
{
    double r = distance(b);
    double f = constants::G * m_mass * b.m_mass / (r * r);
    m_force = {
            m_force.x + f * (b.m_pos.x - m_pos.x) / r,
            m_force.y + f * (b.m_pos.y - m_pos.y) / r};
}
void Body::reset_force()
{
    m_force = {0, 0};
}
void Body::update(double delta_t)
{
    Cartesian a = {m_force.x / m_mass, m_force.y / m_mass};
    m_speed = {
            m_speed.x + delta_t * a.x,
            m_speed.y + delta_t * a.y};
    m_pos = {
            m_pos.x + delta_t * m_speed.x,
            m_pos.y + delta_t * m_speed.y};
    reset_force();
}
bool Body::in(const Quadrant & q) const
{
    return q.contains(m_pos);
}
Body Body::plus(const Body & b) const
{
    double new_mass = m_mass + b.m_mass;
    Cartesian new_pos = {
            (m_pos.x * m_mass + b.m_pos.x * b.m_mass) / new_mass,
            (m_pos.y * m_mass + b.m_pos.y * b.m_mass) / new_mass};
    return Body(new_pos, m_speed, m_mass + b.m_mass, "");
}
std::ostream & operator<<(std::ostream & strm, const Body & value)
{
    return strm << value.m_pos.x << ' ' << value.m_pos.y << ' ' << value.m_speed.x << ' ' << value.m_speed.y << ' ' << value.m_mass << ' ' << value.m_name;
}
const Cartesian & Body::get_m_pos() const
{
    return m_pos;
}
const std::string & Body::get_m_name() const
{
    return m_name;
}

PositionTracker::PositionTracker(const std::string & filename)
{
    std::ifstream file(filename);
    if (!file) {
        throw std::ios::failure("Cannot open " + filename);
    }
    file >> m_size;
    double p_x, p_y, v_x, v_y, m;
    std::string name;
    while (file >> p_x >> p_y >> v_x >> v_y >> m >> name) {
        m_bodies.emplace_back(Cartesian{p_x, p_y}, Cartesian{v_x, v_y}, m, name);
    }
}
Track PositionTracker::track(const std::string & body_name, size_t end_time, size_t time_step)
{
    Track ans;
    Body * target_body = nullptr;
    for (std::size_t from = 0; from < m_bodies.size(); ++from) {
        if (m_bodies[from].get_m_name() == body_name) {
            target_body = &m_bodies[from];
            ans.push_back(target_body->get_m_pos());
            break;
        }
    }

    for (size_t time = time_step; time <= end_time; time += time_step) {
        track_impl();

        for (Body & i : m_bodies) {
            i.update(time_step);
        }
        if (target_body == nullptr) {
            throw std::invalid_argument(body_name + " not found");
        }
        ans.push_back(target_body->get_m_pos());
    }
    return ans;
}
const std::vector<Body> & PositionTracker::bodies()
{
    return m_bodies;
}
BasicPositionTracker::BasicPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}
void BasicPositionTracker::track_impl()
{
    for (Body & from : m_bodies) {
        for (const Body & to : m_bodies) {
            if (from.distance(to) > 0) {
                from.add_force(to);
            }
        }
    }
}

FastPositionTracker::FastPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}
void FastPositionTracker::track_impl()
{
    std::unique_ptr<BHTreeInternal> m_root(new BHTreeInternal({{0, 0}, m_size}));
    for (const auto & i : m_bodies) {
        m_root->insert(i);
    }
    for (Body & from : m_bodies) {
        m_root->update_force(from);
    }
}

BHTreeNode::BHTreeNode(const Body & b)
    : m_representative(b)
{
}
void BHTreeNode::update_force(Body & b)
{
    if (is_far_enough(b)) {
        b.add_force(m_representative);
    }
}
bool BHTreeNode::is_far_enough(const Body & b)
{
    return m_representative.distance(b) > 0;
}
const Body & BHTreeNode::get_m_representative() const
{
    return m_representative;
}
void BHTreeNode::insert(const Body &)
{
    throw std::bad_function_call{};
}

BHTreeInternal::BHTreeInternal(const Quadrant & mArea)
    : m_area(mArea)
{
}
void BHTreeInternal::insert(const Body & b)
{
    m_representative = m_representative.plus(b);

    assert(b.in(m_area));
    std::pair<std::size_t, std::size_t> target;
    Quadrant target_subquadrant = m_area;
    if (b.in(m_area.nw())) {
        target = {0, 0};
        target_subquadrant = m_area.nw();
    }
    else if (b.in(m_area.ne())) {
        target = {1, 0};
        target_subquadrant = m_area.ne();
    }
    else if (b.in(m_area.sw())) {
        target = {0, 1};
        target_subquadrant = m_area.sw();
    }
    else if (b.in(m_area.se())) {
        target = {1, 1};
        target_subquadrant = m_area.se();
    }
    else {
        throw std::logic_error("Body is not in any subquadrant");
    }

    if (m_children[target.first][target.second] == nullptr) {
        m_children[target.first][target.second] = std::make_unique<BHTreeNode>(b);
    }
    else {
        BHTreeInternal * target_internal = dynamic_cast<BHTreeInternal *>(m_children[target.first][target.second].get());
        if (target_internal != nullptr) {
            target_internal->insert(b);
        }
        else {
            Body old_body = m_children[target.first][target.second]->get_m_representative();
            m_children[target.first][target.second] = std::make_unique<BHTreeInternal>(target_subquadrant);
            m_children[target.first][target.second]->insert(old_body);
            m_children[target.first][target.second]->insert(b);
        }
    }
}
bool BHTreeInternal::is_far_enough(const Body & b)
{
    return m_area.length() / b.distance(m_representative) < constants::THETA;
}
void BHTreeInternal::update_force(Body & b)
{
    if (is_far_enough(b)) {
        b.add_force(m_representative);
    }
    else {
        for (const auto & i : m_children) {
            for (const auto & j : i) {
                if (j) {
                    j->update_force(b);
                }
            }
        }
    }
}

Quadrant::Quadrant(Cartesian center, double length)
    : m_center(center)
    , m_length(length)
{
}
bool Quadrant::contains(const Cartesian & p) const
{
    return m_center.x - length() <= p.x && p.x <= m_center.x + length() && m_center.y - length() <= p.y && p.y <= m_center.y + length();
}
double Quadrant::length() const
{
    return m_length;
}
Quadrant Quadrant::nw() const
{
    return Quadrant({m_center.x - length() / 2, m_center.y - length() / 2}, length() / 2);
}
Quadrant Quadrant::ne() const
{
    return Quadrant({m_center.x + length() / 2, m_center.y - length() / 2}, length() / 2);
}
Quadrant Quadrant::sw() const
{
    return Quadrant({m_center.x - length() / 2, m_center.y + length() / 2}, length() / 2);
}
Quadrant Quadrant::se() const
{
    return Quadrant({m_center.x + length() / 2, m_center.y + length() / 2}, length() / 2);
}
std::ostream & operator<<(std::ostream & strm, const Quadrant & value)
{
    return strm << value.m_center.x << ' ' << value.m_center.y << ' ' << value.length();
}
