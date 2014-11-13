#define NUMDIMS 2
#define MAXNODES 4
#define MINNODES 2

/// Minimal bounding rectangle (n-dimensional)
struct Rect
{
	double m_min[NUMDIMS];                      ///< Min dimensions of bounding box 
    double m_max[NUMDIMS];                      ///< Max dimensions of bounding box 
};

struct Node;

/// May be data or may be another subtree
/// The parents level determines this.
/// If the parents level is 0, then this is data
struct Branch
{
	Rect m_rect;                                  ///< Bounds
    union
    {
		Node* m_child;                              ///< Child node
		int m_data;                            ///< Data Id or Ptr
    };
};
/// Node for each branch level
struct Node
{
    bool IsInternalNode()                         { return (m_level > 0); } // Not a leaf, but a internal node
    bool IsLeaf()                                 { return (m_level == 0); } // A leaf, contains data
    
    int m_count;                                  ///< Count
    int m_level;                                  ///< Leaf is zero, others positive
    Branch m_branch[MAXNODES];                    ///< Branch
};
  
/// A link list of nodes for reinsertion after a delete operation
struct ListNode
{
    ListNode* m_next;                             ///< Next in list
    Node* m_node;                                 ///< Node
};

struct K_Nearest
{
	double dist;
	int ID;
};
