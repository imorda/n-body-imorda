#include <cassert>
#include <cmath>
#include <fstream>
#include <nbody.h>

Body::Body(const Cartesian & pos, const Cartesian & speed, const double mass, const std::string & name)
    : m_pos(pos)
    , m_speed(speed)
    , m_mass(mass)
    , m_name(name)
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
bool Body::in(const Quadrant q) const
{
    return q.contains(m_pos);
}
Body Body::plus(const Body & b)
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
const Cartesian & Body::get_position() const
{
    return m_pos;
}
const std::string & Body::get_name() const
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
Track PositionTracker::track_common(const std::string & body_name, size_t end_time, size_t time_step)
{
    Track ans;
    auto target_body = std::find_if(m_bodies.begin(), m_bodies.end(), [&body_name](const Body & body) { return body_name == body.get_name(); });
    if (target_body == m_bodies.end()) {
        throw std::invalid_argument(body_name + " not found");
    }

    for (size_t time = time_step; time <= end_time; time += time_step) {
        track_impl();

        for (Body & i : m_bodies) {
            i.update(time_step);
        }
        ans.push_back(target_body->get_position());
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
Track BasicPositionTracker::track(const std::string & body_name, size_t end_time, size_t time_step)
{
    return track_common(body_name, end_time, time_step);
}

FastPositionTracker::FastPositionTracker(const std::string & filename)
    : PositionTracker(filename)
{
}
void FastPositionTracker::track_impl()
{
    auto default_quadrant = Quadrant({0, 0}, m_size);
    auto m_root = std::make_unique<BHTreeNode>(default_quadrant);
    for (const auto & i : m_bodies) {
        m_root->insert(i, default_quadrant);
    }
    for (Body & from : m_bodies) {
        m_root->update_force(from);
    }
}
Track FastPositionTracker::track(const std::string & body_name, size_t end_time, size_t time_step)
{
    return track_common(body_name, end_time, time_step);
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
    else {
        if (m_internal_data) {
            if (m_internal_data->nw) {
                m_internal_data->nw->update_force(b);
            }
            if (m_internal_data->ne) {
                m_internal_data->ne->update_force(b);
            }
            if (m_internal_data->sw) {
                m_internal_data->sw->update_force(b);
            }
            if (m_internal_data->se) {
                m_internal_data->se->update_force(b);
            }
        }
    }
}
bool BHTreeNode::is_far_enough(const Body & b)
{
    if (m_internal_data) {
        return m_internal_data->area.length() / b.distance(m_representative) < constants::THETA;
    }
    return m_representative.distance(b) > 0;
}
void BHTreeNode::insert(const Body & b, const Quadrant & quadrant)
{
    if (m_internal_data == nullptr) {
        m_internal_data = std::make_unique<InternalData>(quadrant);
        insert_impl(m_representative);
    }
    m_representative = m_representative.plus(b);
    insert_impl(b);
}
void BHTreeNode::insert_impl(const Body & b)
{
    assert(m_internal_data);
    assert(b.in(m_internal_data->area));

    Quadrant target_subquadrant;
    std::unique_ptr<BHTreeNode> * target_child_node;

    if (b.in(m_internal_data->area.nw())) {
        target_subquadrant = m_internal_data->area.nw();
        target_child_node = &m_internal_data->nw;
    }
    else if (b.in(m_internal_data->area.ne())) {
        target_subquadrant = m_internal_data->area.ne();
        target_child_node = &m_internal_data->ne;
    }
    else if (b.in(m_internal_data->area.sw())) {
        target_subquadrant = m_internal_data->area.sw();
        target_child_node = &m_internal_data->sw;
    }
    else if (b.in(m_internal_data->area.se())) {
        target_subquadrant = m_internal_data->area.se();
        target_child_node = &m_internal_data->se;
    }
    else {
        throw std::logic_error("Body is not in any subquadrant");
    }

    if ((*target_child_node) == nullptr) {
        (*target_child_node) = std::make_unique<BHTreeNode>(b);
    }
    else {
        target_child_node->get()->insert(b, target_subquadrant);
    }
}
BHTreeNode::BHTreeNode(const Quadrant & q)
    : m_internal_data(new InternalData(q))
{
}
Quadrant::Quadrant(Cartesian center, double length)
    : m_center(center)
    , m_length(length)
{
}
bool Quadrant::contains(Cartesian p) const
{
    return m_center.x - length() / 2 <= p.x && p.x <= m_center.x + length() / 2 && m_center.y - length() / 2 <= p.y && p.y <= m_center.y + length() / 2;
}
double Quadrant::length() const
{
    return m_length;
}
Quadrant Quadrant::nw() const
{
    return Quadrant({m_center.x - length() / 4, m_center.y - length() / 4}, length() / 2);
}
Quadrant Quadrant::ne() const
{
    return Quadrant({m_center.x + length() / 4, m_center.y - length() / 4}, length() / 2);
}
Quadrant Quadrant::sw() const
{
    return Quadrant({m_center.x - length() / 4, m_center.y + length() / 4}, length() / 2);
}
Quadrant Quadrant::se() const
{
    return Quadrant({m_center.x + length() / 4, m_center.y + length() / 4}, length() / 2);
}
std::ostream & operator<<(std::ostream & strm, const Quadrant & value)
{
    return strm << value.m_center.x << ' ' << value.m_center.y << ' ' << value.length();
}
BHTreeNode::InternalData::InternalData(const Quadrant & area)
    : area(area)
{
}
