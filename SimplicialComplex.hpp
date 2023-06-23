
// note that this code only works for triangle/tet meshes
// TODO: optimize it from n^2 to nlg(n)
struct Simplex {
    int d;
    Tuple t;
    Simplex(int d, Tuple t) : d(d), t(t) {}
};

class SimplicialComplex
{
private:
    std::vector<std::vector<Tuple>> simplexes;
    const Mesh &m;
public:
    std::vector<std::vector<Tuple>> get_simplexes() const
    {
        return simplexes;
    }
    
    int get_size() const
    {
        int ret = 0;
        for (int d = 0; d < 4; d++)
        {
            ret += simplexes[d].size();
        }
        return ret;
    }

    bool AddSimplex(Simplex s)
    {
        assert(s.d < 4);
        for (auto tmp : simplexes[s.d])
        {
            // TODO: need implement is_same_simplex for tuples
            // is_same_simplex(t1, t2, dim, m) returns if (dim, t1) and (dim, t2) are the same simplex in Mesh m
            // Need it in Mesh class
            if (is_same_simplex(tmp, s.t, s.d, m))
            {
                return false;
            }
        }
        simplexes[s.d].push_back(s.t);
        return true;
    }
    
    void unionComplex(const SimplicialComplex &other)
    {
        for (int d = 0; d < 4; d++)
        {
            for (auto t : other.simplexes[d])
            {
                AddSimplex(Simplex(d, t));
            }
        }
    }

    bool operator==(const SimplicialComplex& other) const
    {
        for (int d = 0; d <= 4; d++)
        {
            if (simplexes[d].size() != other.simplexes[d].size())
            {
                return false;
            }   
            for (auto t1 : simplexes[d])
            {
                bool flag = false;
                for (auto t2 : other.simplexes[d])
                {
                    if (is_same_simplex(t1, t2, d, m))
                    {
                        flag = true;
                        break;
                    }
                }
                if (!flag) return false;
            }
        }
        return true;
    }

    SimplicialComplex& operator=(const SimplicialComplex& other) 
    {
        if (this != &other) 
        {
            this->simplexes = other.simplexes;
        }
        return *this;
    }
    
    SimplicialComplex(const Mesh& m);
    ~SimplicialComplex();
};

SimplicialComplex::SimplicialComplex(const Mesh& m)
{
    this->m = m;
    simplexes.resize(4);
}

SimplicialComplex::~SimplicialComplex()
{
}


SimplicialComplex SC_intersect(const SimplicialComplex &A, const SimplicialComplex &B, const Mesh &m)
{
    SimplicialComplex SC_union = A;
    SimplicialComplex SC_intersect(m);
    auto s_B = B.get_simplexes();

    for (int d = 0; d < 4; d++)
    {
        for (auto t : s_B)
        {
            if (!SC_union.AddSimplex(Simplex(d, t)))
            {
                SC_intersect.AddSimplex(Simplex(d, t));
            }
        }
    }

    return SC_intersect;
}

//////////////////////////////////
// List of Operators
// bd: boundary
// clbd: closed boudnary
// st: start
// clst: closed star
// lnk: link
//////////////////////////////////

// ∂s
SimplicialComplex bd(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC(m);
    std::queue<Simplex> q;
    q.push(s);

    while (!q.empty())
    {
        auto cur_s = q.front();
        q.pop();

        switch (cur_s.d - 1)
        {
        case 0: // V
            SC.AddSimplex(Simplex(0, cur_s.t));
            SC.AddSimplex(Simplex(0, cur_s.t.sw(0, m)));
            break;
        
        case 1: // E
            if (SC.AddSimplex(Simplex(1, cur_s.t)))
            {
                q.push(Simplex(1, cur_s.t));
            }
            if (SC.AddSimplex(Simplex(1, cur_s.t.sw(1, m))))
            {
                q.push(Simplex(1, cur_s.t.sw(1, m)));
            }
            if (SC.AddSimplex(Simplex(1, cur_s.t.sw(0, m).sw(1, m))))
            {
                q.push(Simplex(1, cur_s.t.sw(0, m).sw(1, m)));
            }
            break;
            
        case 2: // F
            if (SC.AddSimplex(Simplex(2, cur_s.t)))
            {
                q.push(Simplex(2, cur_s.t));
            }
            if (SC.AddSimplex(Simplex(2, cur_s.t.sw(2, m))))
            {
                q.push(Simplex(2, cur_s.t.sw(2, m)));
            }
            if (SC.AddSimplex(Simplex(2, cur_s.t.sw(1, m).sw(2, m))))
            {
                q.push(Simplex(2, cur_s.t.sw(1, m).sw(2, m)));
            }
            if (SC.AddSimplex(Simplex(2, cur_s.t.sw(0, m).sw(1, m).sw(2, m))))
            {
                q.push(Simplex(2, cur_s.t.sw(0, m).sw(1, m).sw(2, m)));
            }
        }
    }

    return SC;
}

// ∂s∪{s}
SimplicialComplex clbd(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC = bd(s, m);
    SC.AddSimplex(s);
    return SC;
}

// Simplex s1,s2, check if A∩B!=∅
// check is intersect(∂s1, ∂s2) has intersections
bool is_intersect(Simplex s1, Simplex s2, Mesh &m)
{
    SimplicialComplex s1_bd = clbd(s1);
    SimplicialComplex s2_bd = clbd(s2);
    SimplicialComplex s1_s2_int = SC_intersect(s1_bd, s2_bd, m);
    return (s1_s2_int.get_size() != 0);
}

SimplicialComplex clst(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC(m);
    int dim = m.dim; // TODO: 2 for trimesh, 3 for tetmesh need it in Mesh class

    if (dim == 2)
    {
        switch (s.d)
        {
        case 0:
            std::queue<Tuple> q;
            q.push(s.t);
            while (!q.empty())
            {
                auto t = q.front();
                q.pop();
                if (SC.AddSimplex(Simplex(2, t)))
                {
                    if (!t.boundary(m))
                    {
                        q.push(t.sw(2, m));
                    }
                    if (!t.sw(1, m).boundary(m))
                    {
                        q.push(t.sw(1, m).sw(2, m));
                    }
                }
            }
            break;
        case 1:
            SC.addSimplex(Simplex(2, s.t));
            if (!s.t.boundary(m))
            {
                SC.addSimplex(Simplex(2, s.t.sw(2, m)));
            }
            break;
        case 2:
            SC.addSimplex(s);
            break;
        default:
            assert(false);
            break;
        }
        
    }
    else if (dim == 3)
    {
        switch (s.d)
        {
        case 0:
            std::queue<Tuple> q;
            q.push(s.t);
            while (!q.empty())
            {
                auto t = q.front();
                q.pop();
                if (SC.AddSimplex(Simplex(3, t)))
                {
                    if (!t.boundary(m))
                    {
                        q.push(t.sw(3, m));
                    }
                    if (!t.sw(2, m).boundary(m))
                    {
                        q.push(t.sw(2, m).sw(3, m));
                    }
                    if (!t.sw(1, m).sw(2, m).boundary(m))
                    {
                        q.push(t.sw(1, m).sw(2, m).sw(3, m));
                    }
                }
            }
            break;
        case 1:
            std::queue<Tuple> q;
            q.push(s.t);
            while (!q.empty())
            {
                auto t = q.front();
                q.pop();
                if (SC.AddSimplex(Simplex(3, t)))
                {
                    if (!t.boundary(m))
                    {
                        q.push(t.sw(3, m));
                    }
                    if (!t.sw(2, m).boundary(m))
                    {
                        q.push(t.sw(2, m).sw(3, m));
                    }
                }
            }
            break;
        case 2:
            SC.addSimplex(Simplex(3, s.t));
            if (!s.t.boundary(m))
            {
                SC.addSimplex(Simplex(3, s.t.sw(3, m)));
            }
            break;
        case 3:
            SC.addSimplex(s);
            break;
        default:
            assert(false);
            break;
        }  
    }

    auto top_tuples = SC.get_simplexes()[dim];
    for (Tuple t : top_tuples)
    {
      SC.unionComplex(bd(Simplex(dim, t),m));
    }
    return SC;
}

SimplicialComplex lnk(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC_clst = clst(s, m);
    SimplicialComplex SC(m);
    auto simplexes = SC_clst.get_simplexes();
    for (int d = 0; d < 4; d++)
    {
        for (auto t : simplexes[d])
        {
            if (!is_intersect(s, Simplex(d, t)))
            {
                SC.addSimplex(Simplex(d, t));
            }
        }
    }
    return SC;
}

SimplicialComplex st(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC_clst = clst(s, m);
    SimplicialComplex SC;
    SC.addSimplex(s);
    auto simplexes = SC_clst.get_simplexes();
    for (int d = s.d + 1; d < 4; d++)
    {
        for (auto t : simplexes[d])
        {
            if (is_intersect(s, Simplex(d, t)))
            {
                SC.addSimplex(Simplex(d, t));
            }
        }
    }
    return SC;
}


//////////////////////////////////
// check link condition
// input Tuple t --> edge (a,b)
// check if lnk(a) ∩ lnk(b) == lnk(ab)
//////////////////////////////////
bool link_cond(Tuple t, const Mesh& m)
{
    SimplicialComplex lhs = lnk(Simplex(t, 0), m); // lnk(a)
    lhs.unionComplex(lnk(Simplex(t.sw(0, m), 0), m)); // Union lnk(b)

    SimplicialComplex rhs = lnk(Simplex(t, 1), m); // lnk(ab)
    return (lhs == rhs);
}


//////////////////////////////////
// k-ring
//////////////////////////////////
std::vector<Tuple> one_ring(Tuple t, const Mesh& m)
{
    Simplex s(0, t);
    SimplicialComplex s_st = st(s, m);
    std::vector<Tuple> Vs = s_st.get_simplexes()[1];
    for (int i = 0; i < Vs.size(); i++)
    {
        if (is_same_simplex(t, Vs[i], 0, m))
        {
            Vs[i] = Vs[i].sw(0, m);
        }
    }
    return Vs;
}

std::vector<Tuple> k_ring(Tuple t, const Mesh& m, int k)
{
    assert(k >= 1);

    if (k == 1)
    {
        return one_ring(t, m);
    }
    else
    {
        std::vector<Tuple> t_one_ring = one_ring(t, m);
        std::vector<std::vector<Tuple>> all_ts;
        for (auto tmp : t_one_ring)
        {
            all_ts.push_back(k_ring(tmp, m, k - 1));
        }

        std::vector<Tuple> ret;
        for (auto t_vec : all_ts)
        {
            for (auto tmp : t_vec)
            {
                bool flag = true;
                for (int i = 0; i < ret.size(); i++)
                {
                    if (is_same_simplex(ret[i], tmp, 0, m))
                    {
                        flag = false;
                        break;
                    }
                    if (flag)
                    {
                        ret.push_back(tmp);
                    }
                }
            }
        }
    
        return ret;
    }
}


