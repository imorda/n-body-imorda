#include <array>
#include <functional>
#include <iostream>
#include <memory>
#include <string_view>
#include <vector>

namespace constants {
const double G = (6.67e-11);
const double THETA = 0.5;
} // Namespace constants

struct Cartesian
{
    double x = 0;
    double y = 0;
};

// Quadrant representation, required for Problem 2
class Quadrant
{
private:
    Cartesian m_center;
    double m_length = 0;

public:
    // Create quadrant with center (x, y) and size 'length'
    Quadrant() = default;
    Quadrant(Cartesian center, double length);
    // Test if point (x, y) is in the quadrant
    bool contains(Cartesian p) const;

    double length() const;
    // The four methods below construct new Quadrant representing sub-quadrant of the invoking quadrant
    Quadrant nw() const;

    Quadrant ne() const;
    Quadrant sw() const;
    Quadrant se() const;

    friend std::ostream & operator<<(std::ostream &, const Quadrant &);
};

// Single body representation, required for Problem 1 and Problem 2
class Body
{
private:
    Cartesian m_pos;
    Cartesian m_speed;
    Cartesian m_force;
    double m_mass = 0;
    std::string m_name;

public:
    Body(const Cartesian & pos, const Cartesian & speed, double mass, const std::string & name);
    Body() = default;
    double distance(const Body & b) const;
    const Cartesian & get_position() const;
    const std::string & get_name() const;
    // calculate the force-on current body by the 'b' and add the value to accumulated force value
    void add_force(const Body & b);
    // reset accumulated force value
    void reset_force();

    // update body's velocity and position
    void update(double delta_t);

    friend std::ostream & operator<<(std::ostream &, const Body &);

    // The methods below to be done for Burnes-Hut algorithm only
    // Test if body is in quadrant
    bool in(const Quadrant q) const;
    // Create new body representing center-of-m_mass of the invoking body and 'b'
    Body plus(const Body & b);
};

using Track = std::vector<Cartesian>;

class PositionTracker
{
protected:
    double m_size;
    std::vector<Body> m_bodies;
    PositionTracker(const std::string & filename);
    virtual void track_impl() = 0;
    Track track_common(const std::string & body_name, size_t end_time, size_t time_step);

public:
    virtual Track track(const std::string & body_name, size_t end_time, size_t time_step) = 0;
    const std::vector<Body> & bodies();
    virtual ~PositionTracker() = default;
};

class BasicPositionTracker : public PositionTracker
{
protected:
    void track_impl() override;

public:
    BasicPositionTracker(const std::string & filename);
    Track track(const std::string & body_name, size_t end_time, size_t time_step) override;
};

// Burnes-Hut tree representation, required for Problem 2
class BHTreeNode
{
    struct InternalData
    {
        std::unique_ptr<BHTreeNode> ne;
        std::unique_ptr<BHTreeNode> se;
        std::unique_ptr<BHTreeNode> nw;
        std::unique_ptr<BHTreeNode> sw;
        Quadrant area;

        InternalData(const Quadrant & area);
    };

    std::unique_ptr<InternalData> m_internal_data = nullptr;
    Body m_representative;

    bool is_far_enough(const Body & b);
    void insert_impl(const Body &);

public:
    BHTreeNode() = default;
    BHTreeNode(const Body & b);
    BHTreeNode(const Quadrant & q);
    void insert(const Body &, const Quadrant &);
    // Update net acting force-on 'b'
    void update_force(Body & b);
};

class FastPositionTracker : public PositionTracker
{
protected:
    void track_impl() override;

public:
    FastPositionTracker(const std::string & filename);
    Track track(const std::string & body_name, size_t end_time, size_t time_step) override;
};
