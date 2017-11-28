#include <stdlib.h>
#include <stdio.h>

typedef struct {
	double x;
	double y;
} Point;

typedef struct YArrayEntry {
	Point p;
	struct YArrayEntry *left;  // points to y array of the left child
	struct YArrayEntry *right; // points to y array of the right child
} YArrayEntry;

typedef enum {
	INTERNAL_NODE,
	LEAF_NODE,
} NodeType;

typedef struct TreeNode TreeNode;

typedef struct {
	double left_max;
	TreeNode *left;
	TreeNode *right;
	YArrayEntry *y_array;
} InternalNode;

typedef struct {
	Point p;
} LeafNode;

typedef union {
	InternalNode internal_node;
	LeafNode leaf_node;
} NodeUnion;

struct TreeNode {
	NodeType node_type;
	NodeUnion node_union;
} ;

typedef struct {
	Point *point_array;
	size_t size;
} QueryResult;

Point *read_data(char *filename) {
	return NULL;
}

TreeNode *construct_x_tree(Point *input_data) {
	return NULL;
}

QueryResult query(TreeNode *root, Point upper_left, Point upper_right) {
	QueryResult dummy = { NULL, 0};
	return dummy;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		printf("Usage: %s datafile\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	Point *input_data = read_data(argv[1]);

	TreeNode *x_tree = construct_x_tree(input_data);
	free(input_data);

	do {
		Point upper_left;
		Point lower_right;

		printf("Please enter the query (in the order of upper-left.x upper-left.y lower-right.x lower-right.y, all zeros to terminate): ");
		scanf("%lf %lf %lf %lf", &upper_left.x, &upper_left.y, &lower_right.x, &lower_right.y);

		if (upper_left.x == 0 && upper_left.y == 0 & lower_right.x == 0 && lower_right.y == 0)
			break;

		QueryResult query_result = query(x_tree, upper_left, lower_right);

		printf("Query result:\n");
		for (size_t i = 0; i < query_result.size; ++i) {
			printf("%f, %f\n", query_result.point_array[i].x, query_result.point_array[i].x);
		}
		free(query_result.point_array);
	} while (1);

	return 0;
}