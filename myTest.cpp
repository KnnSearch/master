#include "myRTree.h"
#include "TimeCounter.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <vector>
using namespace std;

#define DB2_SIZE MAXNODES 
#define DB3_SIZE 3
#define MAP_SIZE 20
#define DATA_NUMBER 100
#define K_NEAREST 2

struct Line;
struct Leaf;

struct U1
{
	Rect rect;
	int number;
	Leaf *ptr;
};

struct Table1
{
	vector<U1> leaf_node;
};

struct U2
{
	int tag;
	double x;
	double y;
	Line *ptr;	
};

struct Leaf
{
	vector<U2> entry;
};

struct Table2
{
	vector<Leaf> leaf;
};

struct U3
{
	int tag;
	int tail;
};

struct Line
{
	vector<U3> entry;
};

struct Table3
{
	vector<Line> line;
};

//Is element in vector
bool IsContained(vector<int> tag,int k)
{
	for(size_t i = 0; i < tag.size();++i){
		if(tag[i] == k) return true;
	}
	return false;
}

//Is entry in node
bool IsContained(Node *n,vector<int> tag)
{
	for(int i = 0; i < n->m_count; ++i){
		for(size_t j = 0; j < tag.size(); ++j)
			if(n->m_branch[i].m_data == tag[j]) return true;
	}
	return false;
}

//Is node in node_queue
bool IsContained(RTree *tree,Node *node,vector<Node *> node_queue)
{
	Rect r1 = tree->getMBR(node);
	for(size_t i = 0; i < node_queue.size(); ++i){
		Rect r2 = tree->getMBR(node_queue[i]);
		if((r1.m_max[0] == r2.m_max[0])&&(r1.m_max[1] == r2.m_max[1])&&(r1.m_min[0] == r2.m_min[0])&&(r1.m_min[1] == r2.m_min[1]))
			return true;
	}
	return false;
}

//Circle and Rectangle
bool IsCover(Rect r,double q[],double radio)
{
	if(minDist(r,q) > radio) return 0;
	else return 1;
}

double maxDistToNodes(RTree *tree, vector<Node *> n, double q[])
{
	double maxDist = 0.0;
	double temp = 0.0;
	for(size_t i = 0; i < n.size(); ++i){
		Rect r = tree->getMBR(n[i]);
		double card[] = {r.m_min[0],r.m_min[1]}; //one
		temp = objectDist(q,card);
		maxDist = (maxDist > temp)? maxDist : temp;
		card[0] = r.m_max[0]; //two
		temp = objectDist(q,card);
		maxDist = (maxDist > temp)? maxDist : temp;
		card[1] = r.m_max[1]; //three
		temp = objectDist(q,card);
		maxDist = (maxDist > temp)? maxDist : temp;
		card[0] = r.m_min[0]; //four
		temp = objectDist(q,card);
		maxDist = (maxDist > temp)? maxDist : temp;
	}
	return maxDist;
}

int QueryPlan(RTree *tree,int k)
{
	int queryPlan = 0;
	vector<Node *> nodes;
	nodes = (*tree).getNodes(0); //get all leaf nodes
	cout << nodes.size() << endl;

	//calculate for every C
	for(int i = 0; i < MAP_SIZE; ++i){
		for(int j = 0; j < MAP_SIZE; ++j){
//	{{int i = 1,j = 18;
			//STEP 1: retrieve knn entries in a search area C while is an 1*1 square
			vector<int> tag; //retrieved entry tag
			vector<K_Nearest> result;
			for(int m = 0 ;m < 20; m++){
				double d1 = (double)rand()/RAND_MAX;
				double d2 = (double)rand()/RAND_MAX;
				double q[] = {i+d1,j+d2}; //randomly generate 20 points in a square
				tree->kNearestNeighbourSearch(tree->getRoot(),q,&result,k);
				for(size_t i = 0; i < result.size(); ++i){
					if(!IsContained(tag,result[i].ID)){
						tag.push_back(result[i].ID);
					}
				}
				vector<K_Nearest>().swap(result);
			}
			//test
			//for(size_t x = 0;x < tag.size(); ++x){
			//	cout << tag[x] << " ";
			//}
			//cout << endl;
			//getchar();
			//STEP 2: get nodes covering all the retrieved entries
			vector<Node *> node_queue;

			for(size_t s = 0; s < nodes.size();++s){
				if(IsContained(nodes[s],tag))
					node_queue.push_back(nodes[s]);
				
			}

			//test
			//for(size_t x = 0; x < nodes.size(); ++x){
			//	tree->printRec(tree->getMBR(nodes[x]));
			//}
			//cout << node_queue.size() << endl;
			//getchar();
			//STEP 3(1): calculate maxDist
			double maxDist = 0.0;
			{
				double temp = 0.0;
				double vertex1[] = {i,j};
				temp = maxDistToNodes(tree,node_queue,vertex1);
				maxDist = (maxDist > temp) ? maxDist : temp;
				double vertex2[] = {i+1,j};
				temp = maxDistToNodes(tree,node_queue,vertex2);
				maxDist = (maxDist > temp) ? maxDist : temp;
				double vertex3[] = {i,j+1};
				temp = maxDistToNodes(tree,node_queue,vertex3);
				maxDist = (maxDist > temp) ? maxDist : temp;
				double vertex4[] = {i+1,j+1};
				temp = maxDistToNodes(tree,node_queue,vertex4);
				maxDist = (maxDist > temp) ? maxDist : temp;
			}
			//test
			//cout << maxDist << endl;

			//STEP 3(2): Minkowski sum of C, and retrieve overlapped nodes
			//Go around 4 edges
			int count = 0;
			while(count != 5){
				double dx = 0;
				double point[2] = {i+dx,j};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dx = dx + 0.2;
			}
			count = 0;
			while(count != 5){
				double dy = 0;
				double point[2] = {i+1,j+dy};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dy = dy + 0.2;
			}
			count = 0;
			while(count != 5){
				double dy = 0;
				double point[2] = {i,j+dy};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dy = dy + 0.2;
			}
			count = 0;
			while(count != 5){
				double dx = 0;
				double point[2] = {i+dx,j+1};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dx = dx + 0.2;
			}
			//test
			//cout << node_queue.size() << endl;
/*			
			//STEP 4(1): calculate maxDist again
			maxDist = 0.0;
			{
				double temp = 0.0;
				double vertex1[] = {i,j};
				temp = maxDistToNodes(tree,node_queue,vertex1);
				maxDist = (maxDist > temp) ? maxDist : temp;
				double vertex2[] = {i+1,j};
				temp = maxDistToNodes(tree,node_queue,vertex2);
				maxDist = (maxDist > temp) ? maxDist : temp;
				double vertex3[] = {i,j+1};
				temp = maxDistToNodes(tree,node_queue,vertex3);
				maxDist = (maxDist > temp) ? maxDist : temp;
				double vertex4[] = {i+1,j+1};
				temp = maxDistToNodes(tree,node_queue,vertex4);
				maxDist = (maxDist > temp) ? maxDist : temp;
			}
			//test
			//cout << maxDist << endl;

			//STEP 4(2): Minkowski sum of C, and retrieve overlapped nodes again
			//Go around 4 edges
			count = 0;
			while(count != 5){
				double dx = 0;
				double point[2] = {i+dx,j};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dx = dx + 0.2;
			}
			count = 0;
			while(count != 5){
				double dy = 0;
				double point[2] = {i+1,j+dy};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dy = dy + 0.2;
			}
			count = 0;
			while(count != 5){
				double dy = 0;
				double point[2] = {i,j+dy};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dy = dy + 0.2;
			}
			count = 0;
			while(count != 5){
				double dx = 0;
				double point[2] = {i+dx,j+1};
				for(size_t s = 0; s < nodes.size();++s){
					if(!IsContained(tree,nodes[s],node_queue)){ //if not in the node_queue
						if(IsCover(tree->getMBR(nodes[s]),point,maxDist)){
							node_queue.push_back(nodes[s]); //if overlapped,push back
						}
					}
				}
				count++;
				dx = dx + 0.2;
			}
*/			//STEP 5: get max query plan
			queryPlan = ((int)node_queue.size() > queryPlan) ? (int)node_queue.size() : queryPlan;
			//cout <<  i << " " << j << " " << node_queue.size() << endl;

			//clear
			vector<Node *>().swap(node_queue);
			vector<int>().swap(tag);
			vector<K_Nearest>().swap(result);

			//test
			//cout << queryPlan << endl;
			//getchar();
			//exit(1);
		}
	}
	return queryPlan;
}

void generateTable(RTree *tree,Table1 *table1,Table2 *table2, Table3 *table3)
{
	vector<Node *> nodes;
	nodes = (*tree).getNodes(0); //get all leaf nodes

	U1 u1_temp;
	U2 u2_temp;
	U3 u3_temp;
	Leaf leaf_temp;
	Line line_temp;
	int count;

	//generate DB1
	for(size_t i = 0; i < nodes.size(); ++i){
		u1_temp.rect = tree->getMBR(nodes[i]);
		u1_temp.number = nodes[i]->m_count;
		u1_temp.ptr = NULL;
		table1->leaf_node.push_back(u1_temp);
	}
	
	//generate DB2
	for(size_t i = 0; i < nodes.size(); ++i){
		for(int j = 0; j < nodes[i]->m_count; ++j){
			u2_temp.tag = nodes[i]->m_branch[j].m_data;
			u2_temp.x = nodes[i]->m_branch[j].m_rect.m_min[0];
			u2_temp.y = nodes[i]->m_branch[j].m_rect.m_min[1];
			u2_temp.ptr = NULL; //set NULL templately
			leaf_temp.entry.push_back(u2_temp);
		}
		
		//Insert Dummy if necessary
		if(nodes[i]->m_count < DB2_SIZE){
			for(int x = 0; x < DB2_SIZE - nodes[i]->m_count; ++x){
				U2 dummy;
				dummy.ptr = NULL;
				dummy.tag = -1;
				dummy.x = -1;
				dummy.y = -1;
				leaf_temp.entry.push_back(dummy);
			}
		}
		table2->leaf.push_back(leaf_temp);
		vector<U2>().swap(leaf_temp.entry);

		//set DB1 ptr
		table1->leaf_node[i].ptr = &(table2->leaf[i]);
	}

	//set DB1 ptr
	for(size_t i = 0; i < table1->leaf_node.size(); ++i){
		table1->leaf_node[i].ptr = &(table2->leaf[i]);
		//cout << "#" << table1->leaf_node[i].ptr << endl;
		//cout << "#" << &(table2->leaf[i]) << endl;
	}
	
	//generate DB3
	vector<U3> temp;
	for(size_t i = 0; i < table2->leaf.size(); ++i){
		for(int j = 0; j < DB2_SIZE; ++j){
			if(table2->leaf[i].entry[j].tag >= 0){
				u3_temp.tag = table2->leaf[i].entry[j].tag;
				u3_temp.tail = 0;  //tail unset
				temp.push_back(u3_temp);
			}
		}
	}
	int line_number;
	vector<U3>::iterator it;
	it = temp.begin();
	if(temp.size() % DB3_SIZE == 0) line_number = temp.size() / DB3_SIZE;
	else line_number = temp.size() / DB3_SIZE + 1;

	count = temp.size();
	for(int i = 0; i < line_number; ++i){
		for(int j = 0; j < DB3_SIZE; ++j){
			if(count != 0){
				line_temp.entry.push_back(*it);
				it++;
				count--;
			}else{
				U3 dummy;
				dummy.tag = -1;
				dummy.tail = -1;
				line_temp.entry.push_back(dummy);
			}
		}
		table3->line.push_back(line_temp);
		vector<U3>().swap(line_temp.entry);
	}
	//set DB2 ptr
	count = 0;
	for(size_t i = 0; i < table2->leaf.size(); ++i){
		for(int j = 0; j < DB2_SIZE; ++j){
			if(table2->leaf[i].entry[j].tag >= 0){
				table2->leaf[i].entry[j].ptr = &(table3->line[count/DB3_SIZE]);
				count++;
			}
		}
	}
}

void printTable1(Table1 *table1)
{
	cout << "========DATA BASE 1========" << endl;
	cout << left << setw(5) << "NO." << left << setw(25) << "MBR"
		 << left << setw(10) << "Number" << left << setw(10) << "ptr" << endl;
	for(size_t i = 0; i < table1->leaf_node.size(); ++i){
		char temp[50];
		sprintf_s(temp,"(%.2f,%.2f,%.2f,%.2f)",
				table1->leaf_node[i].rect.m_min[0],table1->leaf_node[i].rect.m_min[1],
				table1->leaf_node[i].rect.m_max[0],table1->leaf_node[i].rect.m_max[1]);
		//cout << temp << endl;
		cout << left << setw(5) << i
			 << left << setw(25) << temp 
			 << left << setw(10) << table1->leaf_node[i].number
			 << left << setw(10) << table1->leaf_node[i].ptr
			 << endl;
	}
}

void printTable2(Table2 *table2)
{
	cout << "========DATA BASE 2========" << endl;
	cout << left << setw(5) << "NO." << "ID" << endl;
	for(size_t i = 0; i < table2->leaf.size(); ++i){
		cout << left << setw(5) << i;
		for(int j = 0; j < DB2_SIZE; ++j){
			cout << left << setw(3) << table2->leaf[i].entry[j].tag ;
		}
		cout << endl;
	}
}

void printTable3(Table3 *table3)
{
	cout << "========DATA BASE 3========" << endl;
	cout << left << setw(5) << "NO." << "ID" << endl;
	for(size_t i = 0; i < table3->line.size(); ++i){
		cout << left << setw(5) << i;
		for(int j = 0; j < DB3_SIZE; ++j){
			cout << left << setw(3) << table3->line[i].entry[j].tag;
		}
		cout << endl;
	}
}

void outFileTable1(Table1 *table1)
{
	fstream outfile("database_1.txt",ios::out);
	if(!outfile) exit(1);

	outfile << left << setw(5) << "NO." << left << setw(25) << "MBR"
		 << left << setw(10) << "Number" << left << setw(10) << "ptr" << endl;
	for(size_t i = 0; i < table1->leaf_node.size(); ++i){
		char temp[50];
		sprintf_s(temp,"%.2f %.2f %.2f %.2f",
				table1->leaf_node[i].rect.m_min[0],table1->leaf_node[i].rect.m_min[1],
				table1->leaf_node[i].rect.m_max[0],table1->leaf_node[i].rect.m_max[1]);
		//cout << temp << endl;
		outfile << left << setw(5) << i
			 << left << setw(25) << temp 
			 << left << setw(10) << table1->leaf_node[i].number
			 << left << setw(10) << table1->leaf_node[i].ptr
			 << endl;
	}
}

void outFileTable2(Table2 *table2)
{
	fstream outfile("database_2.txt",ios::out);
	if(!outfile) exit(1);

	outfile << left << setw(5) << "NO." << "ID" << endl;
	for(size_t i = 0; i < table2->leaf.size(); ++i){
		outfile << left << setw(5) << i;
		for(int j = 0; j < DB2_SIZE; ++j){
			outfile << left << setw(3) << table2->leaf[i].entry[j].tag ;
		}
		outfile << endl;
	}
}

void outFileTable3(Table3 *table3)
{
	fstream outfile("database_3.txt",ios::out);
	if(!outfile) exit(1);

	outfile << left << setw(5) << "NO." << "ID" << endl;
	for(size_t i = 0; i < table3->line.size(); ++i){
		outfile << left << setw(5) << i;
		for(int j = 0; j < DB3_SIZE; ++j){
			outfile << left << setw(3) << table3->line[i].entry[j].tag;
		}
		outfile << endl;
	}
}

void printTable(Table1 *table1,Table2 *table2,Table3 *table3)
{
	printTable1(table1);
	printTable2(table2);
	printTable3(table3);
}

void outFileTable(Table1 *table1,Table2 *table2,Table3 *table3)
{
	outFileTable1(table1);
	outFileTable2(table2);
	outFileTable3(table3);
}

void outFileData(RTree *tree)
{
	fstream outfile("data.txt",ios::out);
	if(!outfile) exit(1);

	RTree::Iterator it;
	for( tree->GetFirst(it); 
       !tree->IsNull(it);
       tree->GetNext(it) )
	{
		int value = tree->GetAt(it);
    
		double x,y;
		it.GetCard(&x,&y);
		outfile << value << " " << x << " " << y << endl;
	}
}

void outFileRect(RTree *tree)
{
	fstream outfile("rect.txt",ios::out);
	if(!outfile) exit(1);
	Node *root = tree->getRoot();
	int level = root->m_level;

	while(level >= 0){
		vector<Rect> t;
		t = tree->getRectangles(level);
		for(size_t i = 0; i < t.size(); ++i){
			outfile << t[i].m_min[0] << " " << t[i].m_min[1] << " "
					<< t[i].m_max[0] << " " << t[i].m_max[1] << endl;
		}
		level--;
	}
}


int main()
{
	RTree tree;
	
	for(int i = 0;i < DATA_NUMBER; ++i)
	{
		double min[2],max[2];
		min[0] = rand() % MAP_SIZE;
		max[0] = min[0];
		min[1] = rand() % MAP_SIZE;
		max[1] = min[1];
		tree.Insert(min, max, i);  //Insert (x,y,id)
	}
	outFileData(&tree);
	//tree.fileInsert("data.txt");
	//tree.printRec(tree.getMBR(tree.getRoot()));

	RTree::Iterator it;
	for( tree.GetFirst(it); 
       !tree.IsNull(it);
       tree.GetNext(it) )
	{
		int value = tree.GetAt(it);
    
		double x,y;
		it.GetCard(&x,&y);
		cout << "ID " << value << " : " << "(" << x << "," << y << ")" << endl;
	}

	Table1 table1;
	Table2 table2;
	Table3 table3;
	generateTable(&tree,&table1,&table2,&table3);
	printTable(&table1,&table2,&table3);
	outFileTable(&table1,&table2,&table3);

	TimeCounter tc;
	double CreateTime = 0.0;
	tc.StartTime();
	cout << QueryPlan(&tree,K_NEAREST) << endl;
	CreateTime = tc.EndTime();
	cout << "QueryPlan Time: " << CreateTime << " ms" << endl;

//	outFileRect(&tree);

	system("pause");
	return 0;
}