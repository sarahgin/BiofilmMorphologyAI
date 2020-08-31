from ComputationalBiology.biology_utils.Gene import Gene

def get_gene_id(gene: Gene):
    """

    :type gene: Gene object
    """
    return gene.get_id()


def get_gene_name(gene: Gene):
    return gene.name


def get_type(gene: Gene):
    return gene.type


def get_product_type(gene: Gene):
    return gene.gene_product.type if gene.gene_product is not None else None


def get_strand(gene: Gene):
    return gene.strand


def get_product_description(gene: Gene):
    return gene.gene_product.description if gene.gene_product is not None else None
