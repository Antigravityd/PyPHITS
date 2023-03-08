
class BijectiveTransformer(Transformer):
    """A Lark transformer that can go both ways."""
    def __init__(self, order):
        self.order = order
    def inv_transform(self, out):
