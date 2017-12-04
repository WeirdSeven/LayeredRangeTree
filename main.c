#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <assert.h>

typedef struct {
    double x;
    double y;
} Point;

typedef struct YArrayEntry {
    Point p;                   // the point
    struct YArrayEntry *left;  // points to y array of the left child
    struct YArrayEntry *right; // points to y array of the right child
} YArrayEntry;

typedef struct {
    double left_max;           // max of the left subtree
    double right_max;          // max of the right subtree
    YArrayEntry *y_array;      // y array
    size_t size;               // size of y array
} TreeNode;

typedef struct {
    TreeNode *nodes;
    size_t size;
} Tree;

typedef struct {
    Point *points;
    size_t size;
} PointArray;

void print_point(Point point) {
    double x = point.x;
    double y = point.y;

    if (x == DBL_MAX && y == DBL_MAX) {
        printf("(INFINITY, INFINITY)\n");
    } else if (x == DBL_MAX) {
        printf("(INFINITY, %f)\n", y);
    } else if (y == DBL_MAX) {
        printf("(%f, INFINITY)\n", x);
    } else {
        printf("(%f, %f)\n", x, y);
    }
}

void print_tree_node(TreeNode tree_node) {
    double left_max = tree_node.left_max;
    double right_max = tree_node.right_max;

    if (left_max == DBL_MAX && right_max == DBL_MAX) {
        printf("left_max = [INFINITY], right_max = [INFINITY]\n");
    } else if (left_max == DBL_MAX) {
        assert(0);
        printf("left_max = [INFINITY], right_max = [%f]\n", right_max);
    } else if (right_max == DBL_MAX) {
        printf("left_max = [%f], right_max = [INFINITY]\n", left_max);
    } else {
        printf("left_max = [%f], right_max = [%f]\n", left_max, right_max);
    }

    printf("y array size = [%zu]\n", tree_node.size);
    for (size_t i = 0; i < tree_node.size; ++i) {
        print_point(tree_node.y_array[i].p);
    }
}

void print_tree(Tree tree) {
    for (size_t i = 0; i < tree.size; ++i) {
        printf("Tree node [%zu]: \n", i);
        print_tree_node(tree.nodes[i]);
        printf("\n");
    }
}

void print_point_array(PointArray point_array) {
    for (size_t i = 0; i < point_array.size; ++i) {
        print_point(point_array.points[i]);
    }
}

double pmax(double n1, double n2) {
    if (n1 == DBL_MAX && n2 == DBL_MAX) {
        return DBL_MAX;
    } else if (n1 == DBL_MAX) {
        return n2;
    } else if (n2 == DBL_MAX) {
        return n1;
    } else {
        return n1 > n2 ? n1 : n2;
    }
}

int comparePoint(const void *a, const void *b) {
    Point *point_a = (Point *)a;
    Point *point_b = (Point *)b;
    double a_x = point_a->x;
    double b_x = point_b->x;

    if (a_x < b_x) {
        return -1;
    } else if (a_x == b_x) {
        return 0;
    } else {
        return 1;
    }
}



PointArray read_data(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error: ");
        exit(EXIT_FAILURE);
    }

    size_t num_points = 0;
    fscanf(fp, "%zu", &num_points);

    Point *points = (Point *) malloc(sizeof(Point) * num_points);

    size_t ret;
    double x, y;
    size_t count = 0;
    while ((ret = fscanf(fp, "%lf %lf", &x, &y)) != EOF) {
        points[count].x = x;
        points[count].y = y;
        count += 1;
    }

    fclose(fp);

    if (num_points != count) {
        printf("Error: Wrong number of points.\n");
        exit(EXIT_FAILURE);
    }

    qsort(points, num_points, sizeof(Point), comparePoint);

    PointArray result = { points , num_points };
    return result;
}

Tree construct_x_tree(PointArray input) {
    size_t num_points = input.size;
    size_t height = ceil(log2(num_points));
    size_t num_leaves = exp2(height);
    size_t num_nodes = num_leaves * 2;

    TreeNode *nodes = (TreeNode *) malloc(sizeof(TreeNode) * num_nodes);

    // Constructs real leave nodes
    for (size_t i = 0; i < num_points; ++i) {
        TreeNode *cur_node = &nodes[i + num_leaves];
    
        YArrayEntry *cur_y_array = (YArrayEntry *) malloc(sizeof(YArrayEntry));
        cur_y_array->p = input.points[i];
        cur_y_array->left = cur_y_array->right = NULL;

        cur_node->y_array = cur_y_array;
        cur_node->size = 1;
        cur_node->left_max = cur_node->right_max = input.points[i].x;
    }
    // Constructs dummy leave nodes
    for (size_t i = num_points; i < num_leaves; ++i) {
        TreeNode *cur_node = &nodes[i + num_leaves];

        cur_node->y_array = NULL;
        cur_node->size = 0;
        cur_node->left_max = cur_node->right_max = DBL_MAX;
    }

    // Constructs internal nodes for all levels
    for (size_t i = 1; i <= height; ++i) {
        size_t start = num_leaves >> i;
        size_t end = start << 1;

        // Constructs internal nodes with depth i
        for (size_t j = start; j < end; ++j) {
            TreeNode *cur_node = &nodes[j];
            TreeNode *left_child = &nodes[j << 1];
            TreeNode *right_child = left_child + 1;

            size_t left_size = left_child->size;
            size_t right_size = right_child->size;
            size_t cur_size = left_size + right_size;

            // Construct the y array by merging
            YArrayEntry *cur_y_array;
            // Case I: Both children are empty
            if (cur_size == 0) {
                cur_y_array = NULL;
            } else {
                cur_y_array = (YArrayEntry *) malloc(sizeof(YArrayEntry) * cur_size);

                // Case II: The right child is empty
                if (right_size == 0) {

                    for (size_t k = 0; k < left_size; ++k) {
                        YArrayEntry *cur_entry = &cur_y_array[k];
                        YArrayEntry *left_entry = &left_child->y_array[k];

                        cur_entry->left = left_entry;
                        cur_entry->right = NULL;
                        cur_entry->p = left_entry->p;
                    }
                } // Case III: Both children are non-empty
                else {
                    size_t left_ptr = 0;
                    size_t right_ptr = 0;
                    for (size_t k = 0; k < cur_size; ++k) {
                        YArrayEntry *cur_entry = &cur_y_array[k];

                        if (left_ptr == left_size) {
                            YArrayEntry *left_entry = &left_child->y_array[left_ptr - 1];
                            YArrayEntry *right_entry = &right_child->y_array[right_ptr];

                            cur_entry->left = left_entry;
                            cur_entry->right = right_entry;
                            cur_entry->p = right_entry->p;
                            right_ptr += 1;
                        } else if (right_ptr == right_size) {
                            YArrayEntry *left_entry = &left_child->y_array[left_ptr];
                            YArrayEntry *right_entry = &right_child->y_array[right_ptr - 1];

                            cur_entry->left = left_entry;
                            cur_entry->right = right_entry;
                            cur_entry->p = left_entry->p;
                            left_ptr += 1;
                        } else {
                            YArrayEntry *left_entry = &left_child->y_array[left_ptr];
                            YArrayEntry *right_entry = &right_child->y_array[right_ptr];

                            cur_entry->left = left_entry;
                            cur_entry->right = right_entry;
                            // Sort according to the y value
                            if (left_entry->p.y <= right_entry->p.y) {
                                cur_entry->p = left_entry->p;
                                left_ptr += 1;
                            } else {
                                cur_entry->p = right_entry->p;
                                right_ptr += 1;
                            }
                        }
                    }
                }
            }

            cur_node->y_array = cur_y_array;
            cur_node->size = cur_size;
            cur_node->left_max = pmax(left_child->left_max, left_child->right_max);
            cur_node->right_max = pmax(right_child->left_max, right_child->right_max);
        }
    }

    Tree result = { nodes, num_nodes };
    return result;
}

PointArray query(Tree tree, Point upper_left, Point upper_right) {
    PointArray dummy = { NULL, 0};
    return dummy;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Usage: %s datafile\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    PointArray input = read_data(argv[1]);

    Tree x_tree = construct_x_tree(input);
    free(input.points);

    print_tree(x_tree);

    do {
        Point upper_left;
        Point lower_right;

        printf("Please enter the query (in the order of upper-left.x upper-left.y lower-right.x lower-right.y, all zeros to terminate): ");
        scanf("%lf %lf %lf %lf", &upper_left.x, &upper_left.y, &lower_right.x, &lower_right.y);

        if (upper_left.x == 0 && upper_left.y == 0 & lower_right.x == 0 && lower_right.y == 0)
            break;

        PointArray query_result = query(x_tree, upper_left, lower_right);

        printf("Query result:\n");
        print_point_array(query_result);
        free(query_result.points);
    } while (1);

    return 0;
}