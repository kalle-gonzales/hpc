#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define MAXSIZE 10000

struct node {
	struct node *left;
	struct node *right;
	int value;
	omp_lock_t lock;
};

struct node_allocator {
	struct node *buf;
	struct node *head;
	struct node *end;
	omp_lock_t lock;
};

// Does tree need a lock as well?
struct tree {
	struct node *root;
	omp_lock_t lock;
};

struct node_allocator *new_node_allocator() {
	struct node_allocator *na;
	na = (struct node_allocator *) malloc(sizeof(struct node_allocator));

	na->buf = (struct node *) malloc(sizeof(struct node *) * MAXSIZE);
	na->head = na->buf;
	na->end = &na->buf[MAXSIZE];

	omp_init_lock(&na->lock);

	return na;
}

void init_tree(struct tree *tree) {
	tree->root = NULL;

	omp_init_lock(&tree->lock);
}

// I am not totally sure, if the lock on n is important at the alloc_node. A lock in "insert" should be more important.
// Besides, if we lock na in insert, the lock in this function could be malfunctioning.
struct node *alloc_node(struct node_allocator *na, int value) {
	struct node *n;

	if (na->head >= na->end) {
		return NULL;
	}
	
	n = na->head++;
	n->value = value;
	n->left = NULL;
	n->right = NULL;

	return n;
}

struct node *insert(struct node *node, struct node_allocator *na, int value) {
	
	if (node == NULL) {
		return alloc_node(na, value);
	}
	if (value > node->value) {
		omp_init_lock(&node->lock);
		omp_set_lock(&node->lock);
		node->right = insert(node->right, na, value);
		omp_unset_lock(&node->lock);
	} else if (value < node->value) {
		omp_init_lock(&node->lock);
		omp_set_lock(&node->lock);
		node->left = insert(node->left, na, value);
		omp_unset_lock(&node->lock);
	}
	
	return node;
}

void insert_file(const char *filename, struct tree *tree) {
	FILE *f;
	int value;
	struct node_allocator *na;
	int values[MAXSIZE];
	int count = 0;
	int i;

	f = fopen(filename, "r");
	if (f == NULL) {
		fprintf(stderr, "Could not read file %s.\n", filename);
		return;
	}

	na = new_node_allocator();

	while (count < MAXSIZE && !feof(f) && fscanf(f, "%d ", &value)) {
		values[count++] = value;
	}

	fclose(f);

	/* Wait until all threads have read data-files. */
#pragma omp barrier

	for (i = 0; i < count; ++i) {
		omp_set_lock(&tree->lock);
		tree->root = insert(tree->root, na, values[i]);
		omp_unset_lock(&tree->lock);
	}

	//omp_destroy_lock(&na->lock);
}

int checksum(struct node *node) {
	if (node == NULL) {
		return 0;
	}
	return node->value + checksum(node->left) + checksum(node->right);
}

int main(int argc, char **argv) {
	struct tree tree;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s <input1> [<input2>] [<input3>] ...\n"
				"Each input file creates a new reading thread.\n\n", argv[0]);
		return 1;
	}

	init_tree(&tree);

	omp_set_num_threads(argc - 1);

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		insert_file(argv[tid + 1], &tree);
	}

	printf("checksum: %d\n", checksum(tree.root));

	return 0;
}
