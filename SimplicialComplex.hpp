
// note that this code only works for triangle/tet meshes
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
        assert(s.d <= 4);
        for (auto tmp : simplexes[s.d])
        {
            // TODO: need implement is_same_simplex for tuples
            if (is_same_simplex(tmp, s.t, s.d, m))
            {
                return false;
            }
        }
        simplexes[s.d].push_back(s.t);
        return true;
    }
    
    void union(const SimplicialComplex &other)
    {
        for (int d = 0; d <= 4; d++)
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
                for (auto t2 : other.simplexes[])
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
            this->m = other.m;
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
            if (!SC_union.addSimplex(Simplex(d, t)))
            {
                SC_intersect.addSimplex(Simplex(d, t));
            }
        }
    }

    return SC_intersect;
}


// Simplex s1,s2, check if A∩B!=∅
// check is intersect(∂s1, ∂s2) has intersections
bool is_intersect(Simplex s1, Simplex s2, Mesh &m)
{
    SimplicialComplex s1_bd = clbd(s1);
    SimplicialComplex s2_bd = cldb(s2);
    SimplicialComplex s1_s2_int = SC_intersect(s1_bd, s2_bd, m);
    return (s1_s2_int.get_size() != 0);
}

//////////////////////////////////
// List of Operators
// bd: boundary
// clbd: closed boudnary
// st: start
// clst: closed star
// lnk: link
/////////////////////////////////

// ∂s
SimplicialComplex bd(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC(m);
    for (int d = 0; d < s.d; d++)
    {
        switch (d)
        {
        case 0: // V
            SC.AddSimplex(Simplex(0, s.t));
            SC.AddSimplex(Simplex(0, s.t.sw(0, m)));
            break;
        
        case 1: // E
            SC.AddSimplex(Simplex(1, s.t));
            SC.AddSimplex(Simplex(1, s.t.sw(1, m)));
            SC.AddSimplex(Simplex(1, s.t.sw(0, m).sw(1, m)));
            break;
            
        case 2: // F
            SC.AddSimplex(Simplex(2, s.t));
            SC.AddSimplex(Simplex(2, s.t.sw(2, m)));
            SC.AddSimplex(Simplex(2, s.t.sw(1, m).sw(2, m)));
            SC.AddSimplex(Simplex(2, s.t.sw(0, m).sw(1, m).sw(2, m)));

            break;

        default:
            assert(false);
            break;
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

SimplicialComplex clst(const Simplex &s, const Mesh& m)
{
    SimplicialComplex SC(m);
    int dim = m.dim; // TODO: 2 for trimesh, 3 for tetmesh

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

    auto top_tuples = SC.get_simplexes[dim];
    for (Tuple t : top_tuples)
    {
      SC.union(bd(Simplex(dim, t),m));
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
