class TreeNode:
    def __init__(self, value, left=None, right=None):
        self.value = value
        self.left = left
        self.right = right

    def __repr__(self):
        return f"TreeNode({self.value})"


def build_sum_tree(leaves):
    levels = []  # 每一层是一个列表

    # 初始：创建叶子节点
    current_level = [TreeNode(v) for v in leaves]
    levels.append(current_level)

    # 构建每一层，直到只剩根节点
    while len(current_level) > 1:
        next_level = []
        for i in range(0, len(current_level), 2):
            left = current_level[i]
            right = current_level[i + 1]
            parent = TreeNode(left.value + right.value, left, right)
            next_level.append(parent)
        levels.append(next_level)
        current_level = next_level

    return levels  # 返回一个二维列表，每层一个子列表



# 示例
if __name__ == "__main__":
    leaves = [1, 2, 3, 4, 5, 6, 7, 8]
    all_nodes = build_sum_tree(leaves)

    print("所有节点：")
    for i, level in enumerate(all_nodes):
        print("len level ",i, len(level) )
        print(f"第{i}层:", level)
