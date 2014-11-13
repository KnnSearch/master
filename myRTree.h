#ifndef RTREE_H
#define RTREE_H

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <string>
#include <fstream>
#include "Struct.h" 
using namespace std;

#define ASSERT assert // RTree uses ASSERT( condition )
#ifndef Min
  #define Min __min 
#endif //Min
#ifndef Max
  #define Max __max 
#endif //Max

#define RTREE_DONT_USE_MEMPOOLS // This version does not contain a fixed memory allocator, fill in lines with EXAMPLE to implement one.
#define RTREE_USE_SPHERICAL_VOLUME // Better split classification, may be slower on some systems

double minDist(Rect rect,double card[])
{
	double p[] = {card[0],card[1]};

	double s[] = {rect.m_min[0],rect.m_min[1]};
	double t[] = {rect.m_max[0],rect.m_max[1]};

	double dist = 0;
	double r[2];
	for(int i = 0; i < 2; ++i){
		if(p[i] < s[i]) r[i] = s[i];
		else if(p[i] > t[i]) r[i] = t[i];
		else r[i] = p[i];
		dist += (p[i]-r[i])*(p[i]-r[i]);
	}
	return dist;
}
double minMaxDist(Rect rect,double card[])
{
	double p[] = {card[0],card[1]};

	double s[] = {rect.m_min[0],rect.m_min[1]};
	double t[] = {rect.m_max[0],rect.m_max[1]};

	double dist[2];
	double rmk,rMi;
	int i;
	for(int k = 0;k < 2;++k){
		if(p[k] <= (s[k]+t[k])/2) rmk = s[k];
		else rmk = t[k];
		if(k == 0) i = 1;
		if(k == 1) i = 0;

		if(p[i] >= (s[i]+t[i])/2) rMi = s[i];
		else rMi = t[i];
		dist[k] = (p[k]-rmk)*(p[k]-rmk) + (p[i]-rMi)*(p[i]-rMi);
	}
	return min(dist[0],dist[1]);
}
double objectDist(double q[],double card[])
{
	return (q[0]-card[0])*(q[0]-card[0])+(q[1]-card[1])*(q[1]-card[1]);
}

class RTree
{
public:
  RTree();
  virtual ~RTree(); 
  void Insert(const double a_min[NUMDIMS], const double a_max[NUMDIMS], const int& a_dataId); 
  void Remove(const double a_min[NUMDIMS], const double a_max[NUMDIMS], const int& a_dataId); 
  int Search(const double a_min[NUMDIMS], const double a_max[NUMDIMS], bool __cdecl a_resultCallback(int a_data, void* a_context), void* a_context); 
  void RemoveAll();
  int Count();

  class Iterator
  {
  private:  
    enum { MAX_STACK = 32 }; //  Max stack size. Allows almost n^32 where n is number of branches in node    
    struct StackElement
    {
      Node* m_node;
      int m_branchIndex;
    };
    
  public:  
    Iterator()                                    { Init(); }
    ~Iterator()                                   { }    
    /// Is iterator invalid
    bool IsNull()                                 { return (m_tos <= 0); }
    /// Is iterator pointing to valid data
    bool IsNotNull()                              { return (m_tos > 0); }
    /// Access the current data element. Caller must be sure iterator is not NULL first.
    int& operator*()
    {
      ASSERT(IsNotNull());
      StackElement& curTos = m_stack[m_tos - 1];
      return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
    } 
/*
    /// Access the current data element. Caller must be sure iterator is not NULL first.
    const int& operator*() const
    {
      ASSERT(IsNotNull());
      StackElement& curTos = m_stack[m_tos - 1];
      return curTos.m_node->m_branch[curTos.m_branchIndex].m_data;
    } 
*/
    /// Find the next data element
    bool operator++()                             { return FindNextData(); }

    /// Get the bounds for this node
    void GetBounds(double a_min[NUMDIMS], double a_max[NUMDIMS])
    {
      ASSERT(IsNotNull());
      StackElement& curTos = m_stack[m_tos - 1];

      Branch& curBranch = curTos.m_node->m_branch[curTos.m_branchIndex];
      
      for(int index = 0; index < NUMDIMS; ++index)
      {
        a_min[index] = curBranch.m_rect.m_min[index];
        a_max[index] = curBranch.m_rect.m_max[index];
      }
	}
	void GetCard(double *x,double *y)
	{
		double t1[NUMDIMS];
		double t2[NUMDIMS];
		GetBounds(t1,t2);
		*x = t1[0];
		*y = t1[1];
	}

  private:
  
    /// Reset iterator
    void Init()                                   { m_tos = 0; }
    /// Find the next data element in the tree (For internal use only)
    bool FindNextData()
    {
      for(;;)
      {
        if(m_tos <= 0)
        {
          return false;
        }
        StackElement curTos = Pop(); // Copy stack top cause it may change as we use it

        if(curTos.m_node->IsLeaf())
        {
          // Keep walking through data while we can
          if(curTos.m_branchIndex+1 < curTos.m_node->m_count)
          {
            // There is more data, just point to the next one
            Push(curTos.m_node, curTos.m_branchIndex + 1);
            return true;
          }
          // No more data, so it will fall back to previous level
        }
        else
        {
          if(curTos.m_branchIndex+1 < curTos.m_node->m_count)
          {
            // Push sibling on for future tree walk
            // This is the 'fall back' node when we finish with the current level
            Push(curTos.m_node, curTos.m_branchIndex + 1);
          }
          // Since cur node is not a leaf, push first of next level to get deeper into the tree
          Node* nextLevelnode = curTos.m_node->m_branch[curTos.m_branchIndex].m_child;
          Push(nextLevelnode, 0);
          
          // If we pushed on a new leaf, exit as the data is ready at TOS
          if(nextLevelnode->IsLeaf())
          {
            return true;
          }
        }
      }
    }

    /// Push node and branch onto iteration stack (For internal use only)
    void Push(Node* a_node, int a_branchIndex)
    {
      m_stack[m_tos].m_node = a_node;
      m_stack[m_tos].m_branchIndex = a_branchIndex;
      ++m_tos;
      ASSERT(m_tos <= MAX_STACK);
    }
    
    /// Pop element off iteration stack (For internal use only)
    StackElement& Pop()
    {
      ASSERT(m_tos > 0);
      --m_tos;
      return m_stack[m_tos];
    }

    StackElement m_stack[MAX_STACK];              ///< Stack as we are doing iteration instead of recursion
    int m_tos;                                    ///< Top Of Stack index
  
    friend RTree;
  };

  /// Get 'first' for iteration
  void GetFirst(Iterator& a_it)
  {
    a_it.Init();
    Node* first = m_root;
    while(first)
    {
      if(first->IsInternalNode() && first->m_count > 1)
      {
        a_it.Push(first, 1); // Descend sibling branch later
      }
      else if(first->IsLeaf())
      {
        if(first->m_count)
        {
          a_it.Push(first, 0);
        }
        break;
      }
      first = first->m_branch[0].m_child;
    }
  }  

  /// Get Next for iteration
  void GetNext(Iterator& a_it)                    { ++a_it; }

  /// Is iterator NULL, or at end?
  bool IsNull(Iterator& a_it)                     { return a_it.IsNull(); }

  /// Get object at iterator position
  int& GetAt(Iterator& a_it)                 { return *a_it; }

protected:
  /// Variables for finding a split partition
  struct PartitionVars
  {
    int m_partition[MAXNODES+1];
    int m_total;
    int m_minFill;
    int m_taken[MAXNODES+1];
    int m_count[2];
    Rect m_cover[2];
    double m_area[2];

    Branch m_branchBuf[MAXNODES+1];
    int m_branchCount;
    Rect m_coverSplit;
    double m_coverSplitArea;
  }; 
 
  Node* AllocNode();
  void FreeNode(Node* a_node);
  void InitNode(Node* a_node);
  void InitRect(Rect* a_rect);
  bool InsertRectRec(Rect* a_rect, const int& a_id, Node* a_node, Node** a_newNode, int a_level);
  bool InsertRect(Rect* a_rect, const int& a_id, Node** a_root, int a_level);
  Rect NodeCover(Node* a_node);
  bool AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode);
  void DisconnectBranch(Node* a_node, int a_index);
  int PickBranch(Rect* a_rect, Node* a_node);
  Rect CombineRect(Rect* a_rectA, Rect* a_rectB);
  void SplitNode(Node* a_node, Branch* a_branch, Node** a_newNode);
  double RectSphericalVolume(Rect* a_rect);
  double RectVolume(Rect* a_rect);
  double CalcRectVolume(Rect* a_rect);
  void GetBranches(Node* a_node, Branch* a_branch, PartitionVars* a_parVars);
  void ChoosePartition(PartitionVars* a_parVars, int a_minFill);
  void LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars);
  void InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill);
  void PickSeeds(PartitionVars* a_parVars);
  void Classify(int a_index, int a_group, PartitionVars* a_parVars);
  bool RemoveRect(Rect* a_rect, const int& a_id, Node** a_root);
  bool RemoveRectRec(Rect* a_rect, const int& a_id, Node* a_node, ListNode** a_listNode);
  ListNode* AllocListNode();
  void FreeListNode(ListNode* a_listNode);
  bool Overlap(Rect* a_rectA, Rect* a_rectB);
  void ReInsert(Node* a_node, ListNode** a_listNode);
  bool Search(Node* a_node, Rect* a_rect, int& a_foundCount, bool __cdecl a_resultCallback(int a_data, void* a_context), void* a_context);
  void RemoveAllRec(Node* a_node);
  void Reset();
  void CountRec(Node* a_node, int& a_count);
  
  Node* m_root;                                    ///< Root of tree
  double m_unitSphereVolume;                 ///< Unit sphere constant for required number of dimensions

  ///////////////////////////////////////////-----------------KNN SEARCH--------------------/////////////////////////////////////////////////////////
public:
	void fileInsert(string filename)
	{
		ifstream infile(filename,ios::in);
		if(!infile){
			cerr << "File could not be opened." << endl;
			exit(1);
		}
		double x,y;
		int id;
		while(infile >> id >> x >> y){
			double min[2] = {x,y};
			double max[2] = {x,y};
			Insert(min, max, id);
		}
	}
	Node *getRoot() {return m_root;}
	Rect getMBR(Node *n)
	{
		Rect r = NodeCover(n);
		return r;
	}

	//return all the nodes in level n
	vector<Node *> getNodes(int level) 
	{
		vector<Node *> temp;
		queue<Node *> queue_node;
		Node *node = getRoot();
		queue_node.push(node);
		while(queue_node.front()->m_level != level){
			for(int i = 0; i < queue_node.front()->m_count;++i){
				queue_node.push(queue_node.front()->m_branch[i].m_child);
			}
			queue_node.pop();
		}
		while(queue_node.size() != 0){
			temp.push_back(queue_node.front());
			queue_node.pop();
		}
		return temp;
	} 
	
	//return all the nodes in level n,leaf nodes n = 0
	vector<Rect> getRectangles(int level) 
	{
		vector<Rect> temp;
		vector<Node *> node_temp = getNodes(level);
		for(int i = 0; i < static_cast<int>(node_temp.size()); ++i){
			temp.push_back(getMBR(node_temp[i]));
		}
		return temp;
	}

	void printRec(Rect t)
	{
		cout << t.m_min[0] << " ";
		cout << t.m_min[1] << " ";
		cout << t.m_max[0] << " ";
		cout << t.m_max[1] << endl;
	}
	void kNearestNeighbourSearch(Node *n,double queryPoint[],vector<K_Nearest> *dist_list,int k)
	{
		double q[2] = {queryPoint[0],queryPoint[1]};
		Node *newNode = AllocNode();
		double d;

		vector<Branch> brenchList; //active brench list;

		if(n->IsLeaf())
		{
			for(int i = 0; i < n->m_count; ++i){
				Rect rect = n->m_branch[i].m_rect;
				double r[2] = {rect.m_min[0],rect.m_min[1]};
				d = objectDist(r,q);
				if(static_cast<int>(dist_list->size()) < k)
				{
					K_Nearest temp;
					temp.dist = d;
					temp.ID = n->m_branch[i].m_data;
					dist_list->push_back(temp);
				}
				else
				{
					K_Nearest temp;
					temp.dist = d;
					temp.ID = n->m_branch[i].m_data;
					dist_list->push_back(temp);
					//sort				
					for(int t = 0; t < static_cast<int>(dist_list->size()-1); ++t)
					{
						for(int s = 0; s < static_cast<int>(dist_list->size()-1-t); ++s)
						{
							if((*dist_list)[s].dist > (*dist_list)[s+1].dist)
							{
								temp = (*dist_list)[s];
								(*dist_list)[s] = (*dist_list)[s+1];
								(*dist_list)[s+1] = temp;
							}
						}
					}					
					dist_list->pop_back();
				}
			}	
		}
		else
		{
			//generate a sorted brench list
			genBrenchList(n,queryPoint,&brenchList);

			//test
			//cout << brenchList.size() << endl;
			//for(size_t t = 0; t < brenchList.size(); ++t){
			//	printRec(brenchList[t].m_rect);	
			//}
			//getchar();
			
			//perform downward pruning
			struct group
			{
				double dist;
				int tag;
			};
			vector<group> MM;
			double Kth_MM;
			
			for(int i = 0;i < static_cast<int>(brenchList.size()); ++i){
				group temp;
				temp.dist = minMaxDist(brenchList[i].m_rect,queryPoint);
				temp.tag = i;
				MM.push_back(temp);
			}

			for(int t = 0; t < static_cast<int>(MM.size()-1); ++t)
			{
				for(int s = 0; s < static_cast<int>(MM.size()-1-t); ++s)
				{
					if(MM[s].dist > MM[s+1].dist)
					{
						group temp = MM[s];
						MM[s] = MM[s+1];
						MM[s+1] = temp;
					}
				}
			}

			int count = 0;
			for(int i = 0; i < static_cast<int>(MM.size()); ++i)
			{
				count = count + getEntryUnderNode(brenchList[MM[i].tag].m_child);
				if(count >= k){
					Kth_MM = MM[i].dist;
					break;
				}
			}
			//cout << Kth_MM << endl;
			for(int i = 0;i < static_cast<int>(brenchList.size()); ++i){
				if(minDist(brenchList[i].m_rect,queryPoint) > Kth_MM){
					vector<Branch>::iterator it = brenchList.begin()+i;
					brenchList.erase(it);
				}
			}
			
			//test
			//cout << brenchList.size() << endl;
			//for(size_t t = 0; t < brenchList.size(); ++t){
			//	printRec(brenchList[t].m_rect);	
			//}
			//getchar();
			
			// end downward pruning

			//recursive
			for (int i = 0; i < static_cast<int>(brenchList.size()); ++i){
				newNode = brenchList[i].m_child;
				kNearestNeighbourSearch(newNode,queryPoint,dist_list,k);
			}

			//perform upward pruning
			if( static_cast<int>(dist_list->size()) >= k){
				for (int i = 0; i < static_cast<int>(brenchList.size()); ++i){
					if(minDist(brenchList[i].m_rect,queryPoint) > (*dist_list)[k-1].dist){
						vector<Branch>::iterator it = brenchList.begin()+i;
						brenchList.erase(it);
					}
				}
			}
		}
	}	
	int getEntryUnderNode(Node *n){
		int i = 0;
		entryUnderNode(n,&i);
		return i;
	}
protected:
	void genBrenchList(Node *n,double queryPoint[],vector<Branch> *brenchList)
	{
		int length = n->m_count;
		vector<double> dist;

		// push branches in current node in ABL
		for(int i = 0;i < length; ++i){
			brenchList->push_back(n->m_branch[i]);
		}
		// calculate minDist from Q to each MBR in ABL
		for(int i = 0;i < static_cast<int>(brenchList->size()); ++i){
			dist.push_back(minDist((*brenchList)[i].m_rect,queryPoint));
		}
		//sort based on minDist
		for(int t = 0; t < static_cast<int>(brenchList->size()-1); ++t)
		{
			for(int s = 0; s < static_cast<int>(brenchList->size()-1-t); ++s)
			{
				if(dist[s] > dist[s+1])
				{
					Branch temp;
					temp = (*brenchList)[s];
					(*brenchList)[s] = (*brenchList)[s+1];
					(*brenchList)[s+1] = temp;
				}
			}
		}
	}
	void entryUnderNode(Node *n,int *count)
	{
		if(n->IsInternalNode()){
			for(int i = 0; i < n->m_count ; ++i){
				entryUnderNode(n->m_branch[i].m_child,count);
			}
		}else {
			*count = *count + n->m_count;
		}
	}
	
///////////////////////////////////////////////-----------------END--------------------////////////////////////////////////////////////
};

RTree::RTree()
{
  ASSERT(MAXNODES > MINNODES);
  ASSERT(MINNODES > 0);


  // We only support machine word size simple data type eg. integer index or object pointer.
  // Since we are storing as union with non data branch
  ASSERT(sizeof(int) == sizeof(void*) || sizeof(int) == sizeof(int));

  // Precomputed volumes of the unit spheres for the first few dimensions
  const float UNIT_SPHERE_VOLUMES[] = {
    0.000000f, 2.000000f, 3.141593f, // Dimension  0,1,2
    4.188790f, 4.934802f, 5.263789f, // Dimension  3,4,5
    5.167713f, 4.724766f, 4.058712f, // Dimension  6,7,8
    3.298509f, 2.550164f, 1.884104f, // Dimension  9,10,11
    1.335263f, 0.910629f, 0.599265f, // Dimension  12,13,14
    0.381443f, 0.235331f, 0.140981f, // Dimension  15,16,17
    0.082146f, 0.046622f, 0.025807f, // Dimension  18,19,20 
  };

  m_root = AllocNode();
  m_root->m_level = 0;
  m_unitSphereVolume = (double)UNIT_SPHERE_VOLUMES[NUMDIMS];
}

RTree::~RTree()
{
  Reset(); // Free, or reset node memory
}

void RTree::Insert(const double a_min[NUMDIMS], const double a_max[NUMDIMS], const int& a_dataId)
{
#ifdef _DEBUG
  for(int index=0; index<NUMDIMS; ++index)
  {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect;
  
  for(int axis=0; axis<NUMDIMS; ++axis)
  {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }
  
  InsertRect(&rect, a_dataId, &m_root, 0);
}

void RTree::Remove(const double a_min[NUMDIMS], const double a_max[NUMDIMS], const int& a_dataId)
{
#ifdef _DEBUG
  for(int index=0; index<NUMDIMS; ++index)
  {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect;
  
  for(int axis=0; axis<NUMDIMS; ++axis)
  {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }

  RemoveRect(&rect, a_dataId, &m_root);
}

int RTree::Search(const double a_min[NUMDIMS], const double a_max[NUMDIMS], bool __cdecl a_resultCallback(int a_data, void* a_context), void* a_context)
{
#ifdef _DEBUG
  for(int index=0; index<NUMDIMS; ++index)
  {
    ASSERT(a_min[index] <= a_max[index]);
  }
#endif //_DEBUG

  Rect rect;
  
  for(int axis=0; axis<NUMDIMS; ++axis)
  {
    rect.m_min[axis] = a_min[axis];
    rect.m_max[axis] = a_max[axis];
  }

  // NOTE: May want to return search result another way, perhaps returning the number of found elements here.

  int foundCount = 0;
  Search(m_root, &rect, foundCount, a_resultCallback, a_context);

  return foundCount;
}

int RTree::Count()
{
  int count = 0;
  CountRec(m_root, count);
  
  return count;
}

void RTree::CountRec(Node* a_node, int& a_count)
{
  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      CountRec(a_node->m_branch[index].m_child, a_count);
    }
  }
  else // A leaf node
  {
    a_count += a_node->m_count;
  }
}

void RTree::RemoveAll()
{
  // Delete all existing nodes
  Reset();

  m_root = AllocNode();
  m_root->m_level = 0;
}

void RTree::Reset()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
  // Delete all existing nodes
  RemoveAllRec(m_root);
#else // RTREE_DONT_USE_MEMPOOLS
  // Just reset memory pools.  We are not using complex types
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

void RTree::RemoveAllRec(Node* a_node)
{
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);

  if(a_node->IsInternalNode()) // This is an internal node in the tree
  {
    for(int index=0; index < a_node->m_count; ++index)
    {
      RemoveAllRec(a_node->m_branch[index].m_child);
    }
  }
  FreeNode(a_node); 
}

Node* RTree::AllocNode()
{
  Node* newNode;
#ifdef RTREE_DONT_USE_MEMPOOLS
  newNode = new Node;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
  InitNode(newNode);
  return newNode;
}

void RTree::FreeNode(Node* a_node)
{
  ASSERT(a_node);

#ifdef RTREE_DONT_USE_MEMPOOLS
  delete a_node;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

// Allocate space for a node in the list used in DeletRect to
// store Nodes that are too empty.

ListNode* RTree::AllocListNode()
{
#ifdef RTREE_DONT_USE_MEMPOOLS
  return new ListNode;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

void RTree::FreeListNode(ListNode* a_listNode)
{
#ifdef RTREE_DONT_USE_MEMPOOLS
  delete a_listNode;
#else // RTREE_DONT_USE_MEMPOOLS
  // EXAMPLE
#endif // RTREE_DONT_USE_MEMPOOLS
}

void RTree::InitNode(Node* a_node)
{
  a_node->m_count = 0;
  a_node->m_level = -1;
}

void RTree::InitRect(Rect* a_rect)
{
  for(int index = 0; index < NUMDIMS; ++index)
  {
    a_rect->m_min[index] = (double)0;
    a_rect->m_max[index] = (double)0;
  }
}

// Inserts a new data rectangle into the index structure.
// Recursively descends tree, propagates splits back up.
// Returns 0 if node was not split.  Old node updated.
// If node was split, returns 1 and sets the pointer pointed to by
// new_node to point to the new node.  Old node updated to become one of two.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.

bool RTree::InsertRectRec(Rect* a_rect, const int& a_id, Node* a_node, Node** a_newNode, int a_level)
{
  ASSERT(a_rect && a_node && a_newNode);
  ASSERT(a_level >= 0 && a_level <= a_node->m_level);

  int index;
  Branch branch;
  Node* otherNode;

  // Still above level for insertion, go down tree recursively
  if(a_node->m_level > a_level)
  {
    index = PickBranch(a_rect, a_node);
    if (!InsertRectRec(a_rect, a_id, a_node->m_branch[index].m_child, &otherNode, a_level))
    {
      // Child was not split
      a_node->m_branch[index].m_rect = CombineRect(a_rect, &(a_node->m_branch[index].m_rect));
      return false;
    }
    else // Child was split
    {
      a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
      branch.m_child = otherNode;
      branch.m_rect = NodeCover(otherNode);
      return AddBranch(&branch, a_node, a_newNode);
    }
  }
  else if(a_node->m_level == a_level) // Have reached level for insertion. Add rect, split if necessary
  {
    branch.m_rect = *a_rect;
    branch.m_child = (Node*) a_id;
    // Child field of leaves contains id of data record
    return AddBranch(&branch, a_node, a_newNode);
  }
  else
  {
    // Should never occur
    ASSERT(0);
    return false;
  }
}

// Insert a data rectangle into an index structure.
// InsertRect provides for splitting the root;
// returns 1 if root was split, 0 if it was not.
// The level argument specifies the number of steps up from the leaf
// level to insert; e.g. a data rectangle goes in at level = 0.
// InsertRect2 does the recursion.
//

bool RTree::InsertRect(Rect* a_rect, const int& a_id, Node** a_root, int a_level)
{
  ASSERT(a_rect && a_root);
  ASSERT(a_level >= 0 && a_level <= (*a_root)->m_level);
#ifdef _DEBUG
  for(int index=0; index < NUMDIMS; ++index)
  {
    ASSERT(a_rect->m_min[index] <= a_rect->m_max[index]);
  }
#endif //_DEBUG  

  Node* newRoot;
  Node* newNode;
  Branch branch;

  if(InsertRectRec(a_rect, a_id, *a_root, &newNode, a_level))  // Root split
  {
    newRoot = AllocNode();  // Grow tree taller and new root
    newRoot->m_level = (*a_root)->m_level + 1;
    branch.m_rect = NodeCover(*a_root);
    branch.m_child = *a_root;
    AddBranch(&branch, newRoot, NULL);
    branch.m_rect = NodeCover(newNode);
    branch.m_child = newNode;
    AddBranch(&branch, newRoot, NULL);
    *a_root = newRoot;
    return true;
  }

  return false;
}

// Find the smallest rectangle that includes all rectangles in branches of a node.

Rect RTree::NodeCover(Node* a_node)
{
  ASSERT(a_node);
  
  int firstTime = true;
  Rect rect;
  InitRect(&rect);
  
  for(int index = 0; index < a_node->m_count; ++index)
  {
    if(firstTime)
    {
      rect = a_node->m_branch[index].m_rect;
      firstTime = false;
    }
    else
    {
      rect = CombineRect(&rect, &(a_node->m_branch[index].m_rect));
    }
  }
  
  return rect;
}


// Add a branch to a node.  Split the node if necessary.
// Returns 0 if node not split.  Old node updated.
// Returns 1 if node split, sets *new_node to address of new node.
// Old node updated, becomes one of two.

bool RTree::AddBranch(Branch* a_branch, Node* a_node, Node** a_newNode)
{
  ASSERT(a_branch);
  ASSERT(a_node);

  if(a_node->m_count < MAXNODES)  // Split won't be necessary
  {
    a_node->m_branch[a_node->m_count] = *a_branch;
    ++a_node->m_count;

    return false;
  }
  else
  {
    ASSERT(a_newNode);
    
    SplitNode(a_node, a_branch, a_newNode);
    return true;
  }
}


// Disconnect a dependent node.
// Caller must return (or stop using iteration index) after this as count has changed

void RTree::DisconnectBranch(Node* a_node, int a_index)
{
  ASSERT(a_node && (a_index >= 0) && (a_index < MAXNODES));
  ASSERT(a_node->m_count > 0);

  // Remove element by swapping with the last element to prevent gaps in array
  a_node->m_branch[a_index] = a_node->m_branch[a_node->m_count - 1];
  
  --a_node->m_count;
}


// Pick a branch.  Pick the one that will need the smallest increase
// in area to accomodate the new rectangle.  This will result in the
// least total area for the covering rectangles in the current node.
// In case of a tie, pick the one which was smaller before, to get
// the best resolution when searching.
int RTree::PickBranch(Rect* a_rect, Node* a_node)
{
  ASSERT(a_rect && a_node);
  
  bool firstTime = true;
  double increase;
  double bestIncr = (double)-1;
  double area;
  double bestArea;
  int best;
  Rect tempRect;

  for(int index=0; index < a_node->m_count; ++index)
  {
    Rect* curRect = &a_node->m_branch[index].m_rect;
    area = CalcRectVolume(curRect);
    tempRect = CombineRect(a_rect, curRect);
    increase = CalcRectVolume(&tempRect) - area;
    if((increase < bestIncr) || firstTime)
    {
      best = index;
      bestArea = area;
      bestIncr = increase;
      firstTime = false;
    }
    else if((increase == bestIncr) && (area < bestArea))
    {
      best = index;
      bestArea = area;
      bestIncr = increase;
    }
  }
  return best;
}

// Combine two rectangles into larger one containing both
Rect RTree::CombineRect(Rect* a_rectA, Rect* a_rectB)
{
  ASSERT(a_rectA && a_rectB);

  Rect newRect;

  for(int index = 0; index < NUMDIMS; ++index)
  {
    newRect.m_min[index] = Min(a_rectA->m_min[index], a_rectB->m_min[index]);
    newRect.m_max[index] = Max(a_rectA->m_max[index], a_rectB->m_max[index]);
  }

  return newRect;
}

// Split a node.
// Divides the nodes branches and the extra one between two nodes.
// Old node is one of the new ones, and one really new one is created.
// Tries more than one method for choosing a partition, uses best result.
void RTree::SplitNode(Node* a_node, Branch* a_branch, Node** a_newNode)
{
  ASSERT(a_node);
  ASSERT(a_branch);

  // Could just use local here, but member or external is faster since it is reused
  PartitionVars localVars;
  PartitionVars* parVars = &localVars;
  int level;

  // Load all the branches into a buffer, initialize old node
  level = a_node->m_level;
  GetBranches(a_node, a_branch, parVars);

  // Find partition
  ChoosePartition(parVars, MINNODES);

  // Put branches from buffer into 2 nodes according to chosen partition
  *a_newNode = AllocNode();
  (*a_newNode)->m_level = a_node->m_level = level;
  LoadNodes(a_node, *a_newNode, parVars);
  
  ASSERT((a_node->m_count + (*a_newNode)->m_count) == parVars->m_total);
}


// Calculate the n-dimensional volume of a rectangle
double RTree::RectVolume(Rect* a_rect)
{
  ASSERT(a_rect);
  
  double volume = (double)1;

  for(int index=0; index<NUMDIMS; ++index)
  {
    volume *= a_rect->m_max[index] - a_rect->m_min[index];
  }
  
  ASSERT(volume >= (double)0);
  
  return volume;
}

// The exact volume of the bounding sphere for the given Rect
double RTree::RectSphericalVolume(Rect* a_rect)
{
  ASSERT(a_rect);
   
  double sumOfSquares = (double)0;
  double radius;

  for(int index=0; index < NUMDIMS; ++index) 
  {
    double halfExtent = ((double)a_rect->m_max[index] - (double)a_rect->m_min[index]) * 0.5f;
    sumOfSquares += halfExtent * halfExtent;
  }

  radius = (double)sqrt(sumOfSquares);
  
  // Pow maybe slow, so test for common dims like 2,3 and just use x*x, x*x*x.
  if(NUMDIMS == 3)
  {
    return (radius * radius * radius * m_unitSphereVolume);
  }
  else if(NUMDIMS == 2)
  {
    return (radius * radius * m_unitSphereVolume);
  }
  else
  {
    return (double)(pow(radius, NUMDIMS) * m_unitSphereVolume);
  }
}

// Use one of the methods to calculate retangle volume
double RTree::CalcRectVolume(Rect* a_rect)
{
#ifdef RTREE_USE_SPHERICAL_VOLUME
  return RectSphericalVolume(a_rect); // Slower but helps certain merge cases
#else // RTREE_USE_SPHERICAL_VOLUME
  return RectVolume(a_rect); // Faster but can cause poor merges
#endif // RTREE_USE_SPHERICAL_VOLUME  
}

// Load branch buffer with branches from full node plus the extra branch.
void RTree::GetBranches(Node* a_node, Branch* a_branch, PartitionVars* a_parVars)
{
  ASSERT(a_node);
  ASSERT(a_branch);

  ASSERT(a_node->m_count == MAXNODES);
    
  // Load the branch buffer
  for(int index=0; index < MAXNODES; ++index)
  {
    a_parVars->m_branchBuf[index] = a_node->m_branch[index];
  }
  a_parVars->m_branchBuf[MAXNODES] = *a_branch;
  a_parVars->m_branchCount = MAXNODES + 1;

  // Calculate rect containing all in the set
  a_parVars->m_coverSplit = a_parVars->m_branchBuf[0].m_rect;
  for(int index=1; index < MAXNODES+1; ++index)
  {
    a_parVars->m_coverSplit = CombineRect(&a_parVars->m_coverSplit, &a_parVars->m_branchBuf[index].m_rect);
  }
  a_parVars->m_coverSplitArea = CalcRectVolume(&a_parVars->m_coverSplit);

  InitNode(a_node);
}

// Method #0 for choosing a partition:
// As the seeds for the two groups, pick the two rects that would waste the
// most area if covered by a single rectangle, i.e. evidently the worst pair
// to have in the same group.
// Of the remaining, one at a time is chosen to be put in one of the two groups.
// The one chosen is the one with the greatest difference in area expansion
// depending on which group - the rect most strongly attracted to one group
// and repelled from the other.
// If one group gets too full (more would force other group to violate min
// fill requirement) then other group gets the rest.
// These last are the ones that can go in either group most easily.
void RTree::ChoosePartition(PartitionVars* a_parVars, int a_minFill)
{
  ASSERT(a_parVars);
  
  double biggestDiff;
  int group, chosen, betterGroup;
  
  InitParVars(a_parVars, a_parVars->m_branchCount, a_minFill);
  PickSeeds(a_parVars);

  while (((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
       && (a_parVars->m_count[0] < (a_parVars->m_total - a_parVars->m_minFill))
       && (a_parVars->m_count[1] < (a_parVars->m_total - a_parVars->m_minFill)))
  {
    biggestDiff = (double) -1;
    for(int index=0; index<a_parVars->m_total; ++index)
    {
      if(!a_parVars->m_taken[index])
      {
        Rect* curRect = &a_parVars->m_branchBuf[index].m_rect;
        Rect rect0 = CombineRect(curRect, &a_parVars->m_cover[0]);
        Rect rect1 = CombineRect(curRect, &a_parVars->m_cover[1]);
        double growth0 = CalcRectVolume(&rect0) - a_parVars->m_area[0];
        double growth1 = CalcRectVolume(&rect1) - a_parVars->m_area[1];
        double diff = growth1 - growth0;
        if(diff >= 0)
        {
          group = 0;
        }
        else
        {
          group = 1;
          diff = -diff;
        }

        if(diff > biggestDiff)
        {
          biggestDiff = diff;
          chosen = index;
          betterGroup = group;
        }
        else if((diff == biggestDiff) && (a_parVars->m_count[group] < a_parVars->m_count[betterGroup]))
        {
          chosen = index;
          betterGroup = group;
        }
      }
    }
    Classify(chosen, betterGroup, a_parVars);
  }

  // If one group too full, put remaining rects in the other
  if((a_parVars->m_count[0] + a_parVars->m_count[1]) < a_parVars->m_total)
  {
    if(a_parVars->m_count[0] >= a_parVars->m_total - a_parVars->m_minFill)
    {
      group = 1;
    }
    else
    {
      group = 0;
    }
    for(int index=0; index<a_parVars->m_total; ++index)
    {
      if(!a_parVars->m_taken[index])
      {
        Classify(index, group, a_parVars);
      }
    }
  }

  ASSERT((a_parVars->m_count[0] + a_parVars->m_count[1]) == a_parVars->m_total);
  ASSERT((a_parVars->m_count[0] >= a_parVars->m_minFill) && 
        (a_parVars->m_count[1] >= a_parVars->m_minFill));
}


// Copy branches from the buffer into two nodes according to the partition.
void RTree::LoadNodes(Node* a_nodeA, Node* a_nodeB, PartitionVars* a_parVars)
{
  ASSERT(a_nodeA);
  ASSERT(a_nodeB);
  ASSERT(a_parVars);

  for(int index=0; index < a_parVars->m_total; ++index)
  {
    ASSERT(a_parVars->m_partition[index] == 0 || a_parVars->m_partition[index] == 1);
    
    if(a_parVars->m_partition[index] == 0)
    {
      AddBranch(&a_parVars->m_branchBuf[index], a_nodeA, NULL);
    }
    else if(a_parVars->m_partition[index] == 1)
    {
      AddBranch(&a_parVars->m_branchBuf[index], a_nodeB, NULL);
    }
  }
}

// Initialize a PartitionVars structure.
void RTree::InitParVars(PartitionVars* a_parVars, int a_maxRects, int a_minFill)
{
  ASSERT(a_parVars);

  a_parVars->m_count[0] = a_parVars->m_count[1] = 0;
  a_parVars->m_area[0] = a_parVars->m_area[1] = (double)0;
  a_parVars->m_total = a_maxRects;
  a_parVars->m_minFill = a_minFill;
  for(int index=0; index < a_maxRects; ++index)
  {
    a_parVars->m_taken[index] = false;
    a_parVars->m_partition[index] = -1;
  }
}

void RTree::PickSeeds(PartitionVars* a_parVars)
{
  int seed0, seed1;
  double worst, waste;
  double area[MAXNODES+1];

  for(int index=0; index<a_parVars->m_total; ++index)
  {
    area[index] = CalcRectVolume(&a_parVars->m_branchBuf[index].m_rect);
  }

  worst = -a_parVars->m_coverSplitArea - 1;
  for(int indexA=0; indexA < a_parVars->m_total-1; ++indexA)
  {
    for(int indexB = indexA+1; indexB < a_parVars->m_total; ++indexB)
    {
      Rect oneRect = CombineRect(&a_parVars->m_branchBuf[indexA].m_rect, &a_parVars->m_branchBuf[indexB].m_rect);
      waste = CalcRectVolume(&oneRect) - area[indexA] - area[indexB];
      if(waste > worst)
      {
        worst = waste;
        seed0 = indexA;
        seed1 = indexB;
      }
    }
  }
  Classify(seed0, 0, a_parVars);
  Classify(seed1, 1, a_parVars);
}

// Put a branch in one of the groups.
void RTree::Classify(int a_index, int a_group, PartitionVars* a_parVars)
{
  ASSERT(a_parVars);
  ASSERT(!a_parVars->m_taken[a_index]);

  a_parVars->m_partition[a_index] = a_group;
  a_parVars->m_taken[a_index] = true;

  if (a_parVars->m_count[a_group] == 0)
  {
    a_parVars->m_cover[a_group] = a_parVars->m_branchBuf[a_index].m_rect;
  }
  else
  {
    a_parVars->m_cover[a_group] = CombineRect(&a_parVars->m_branchBuf[a_index].m_rect, &a_parVars->m_cover[a_group]);
  }
  a_parVars->m_area[a_group] = CalcRectVolume(&a_parVars->m_cover[a_group]);
  ++a_parVars->m_count[a_group];
}

// Delete a data rectangle from an index structure.
// Pass in a pointer to a Rect, the tid of the record, ptr to ptr to root node.
// Returns 1 if record not found, 0 if success.
// RemoveRect provides for eliminating the root.

bool RTree::RemoveRect(Rect* a_rect, const int& a_id, Node** a_root)
{
  ASSERT(a_rect && a_root);
  ASSERT(*a_root);

  Node* tempNode;
  ListNode* reInsertList = NULL;

  if(!RemoveRectRec(a_rect, a_id, *a_root, &reInsertList))
  {
    // Found and deleted a data item
    // Reinsert any branches from eliminated nodes
    while(reInsertList)
    {
      tempNode = reInsertList->m_node;

      for(int index = 0; index < tempNode->m_count; ++index)
      {
        InsertRect(&(tempNode->m_branch[index].m_rect),
                   tempNode->m_branch[index].m_data,
                   a_root,
                   tempNode->m_level);
      }
      
      ListNode* remLNode = reInsertList;
      reInsertList = reInsertList->m_next;
      
      FreeNode(remLNode->m_node);
      FreeListNode(remLNode);
    }
    
    // Check for redundant root (not leaf, 1 child) and eliminate
    if((*a_root)->m_count == 1 && (*a_root)->IsInternalNode())
    {
      tempNode = (*a_root)->m_branch[0].m_child;
      
      ASSERT(tempNode);
      FreeNode(*a_root);
      *a_root = tempNode;
    }
    return false;
  }
  else
  {
    return true;
  }
}

// Delete a rectangle from non-root part of an index structure.
// Called by RemoveRect.  Descends tree recursively,
// merges branches on the way back up.
// Returns 1 if record not found, 0 if success.

bool RTree::RemoveRectRec(Rect* a_rect, const int& a_id, Node* a_node, ListNode** a_listNode)
{
  ASSERT(a_rect && a_node && a_listNode);
  ASSERT(a_node->m_level >= 0);

  if(a_node->IsInternalNode())  // not a leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      if(Overlap(a_rect, &(a_node->m_branch[index].m_rect)))
      {
        if(!RemoveRectRec(a_rect, a_id, a_node->m_branch[index].m_child, a_listNode))
        {
          if(a_node->m_branch[index].m_child->m_count >= MINNODES)
          {
            // child removed, just resize parent rect
            a_node->m_branch[index].m_rect = NodeCover(a_node->m_branch[index].m_child);
          }
          else
          {
            // child removed, not enough entries in node, eliminate node
            ReInsert(a_node->m_branch[index].m_child, a_listNode);
            DisconnectBranch(a_node, index); // Must return after this call as count has changed
          }
          return false;
        }
      }
    }
    return true;
  }
  else // A leaf node
  {
    for(int index = 0; index < a_node->m_count; ++index)
    {
      if(a_node->m_branch[index].m_child == (Node*)a_id)
      {
        DisconnectBranch(a_node, index); // Must return after this call as count has changed
        return false;
      }
    }
    return true;
  }
}

// Decide whether two rectangles overlap.
bool RTree::Overlap(Rect* a_rectA, Rect* a_rectB)
{
  ASSERT(a_rectA && a_rectB);

  for(int index=0; index < NUMDIMS; ++index)
  {
    if (a_rectA->m_min[index] > a_rectB->m_max[index] ||
        a_rectB->m_min[index] > a_rectA->m_max[index])
    {
      return false;
    }
  }
  return true;
}

// Add a node to the reinsertion list.  All its branches will later
// be reinserted into the index structure.
void RTree::ReInsert(Node* a_node, ListNode** a_listNode)
{
  ListNode* newListNode;

  newListNode = AllocListNode();
  newListNode->m_node = a_node;
  newListNode->m_next = *a_listNode;
  *a_listNode = newListNode;
}

// Search in an index tree or subtree for all data retangles that overlap the argument rectangle.
bool RTree::Search(Node* a_node, Rect* a_rect, int& a_foundCount, bool __cdecl a_resultCallback(int a_data, void* a_context), void* a_context)
{
  ASSERT(a_node);
  ASSERT(a_node->m_level >= 0);
  ASSERT(a_rect);

  if(a_node->IsInternalNode()) // This is an internal node in the tree
  {
    for(int index=0; index < a_node->m_count; ++index)
    {
      if(Overlap(a_rect, &a_node->m_branch[index].m_rect))
      {
        if(!Search(a_node->m_branch[index].m_child, a_rect, a_foundCount, a_resultCallback, a_context))
        {
          return false; // Don't continue searching
        }
      }
    }
  }
  else // This is a leaf node
  {
    for(int index=0; index < a_node->m_count; ++index)
    {
      if(Overlap(a_rect, &a_node->m_branch[index].m_rect))
      {
        int& id = a_node->m_branch[index].m_data;
        
        // NOTE: There are different ways to return results.  Here's where to modify
        if(&a_resultCallback)
        {
          ++a_foundCount;
          if(!a_resultCallback(id, a_context))
          {
            return false; // Don't continue searching
          }
        }
      }
    }
  }

  return true; // Continue searching
}

#undef RTree

#endif //RTREE_H