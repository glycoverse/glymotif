# ========== Glycan Graphs ==========
make_glycan <- function(iupac, mono_type, linkage) {
  glycan <- glyparse::parse_iupac_condensed(iupac)
  glycan <- glyrepr::convert_mono_type(glycan, to = mono_type)
  if (!linkage) {
    glycan <- glyrepr::remove_linkages(glycan)
  }
  glycan
}

paucimannose_H3N2 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     └─Man (a1-3)
  make_glycan("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

paucimannose_H4N2a3 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-3)
  #     └─Man (a1-6)
  #       └─Man (a1-3)
  make_glycan("Man(a1-3)Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

paucimannose_H4N2a6 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-3)
  #     └─Man (a1-6)
  #       └─Man (a1-6)
  make_glycan("Man(a1-6)Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

highmannose_H5N2 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-3)
  #     └─Man (a1-6)
  #       ├─Man (a1-6)
  #       └─Man (a1-3)
  make_glycan("Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

highmannose_H6N2 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─Man (a1-2)
  make_glycan("Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

hybrid_H5N3 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

hybrid_H6N3 <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  make_glycan("Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

hybrid_H4N3a3 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

hybrid_H4N3a6 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-6)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

hybrid_H4N3a3F1 <- function(mono_type, linkage) {
  # GlcNAc
  # ├─Fuc (a1-6)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─Man (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-", mono_type, linkage)
}

complex_H3N3 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H3N4 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H4N4 <- function(mono_type, linkage) {
  # GlcNAc
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  make_glycan("Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H4N4S1 <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  #           └─Neu5Ac (a2-3)
  make_glycan("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H3N5_bisect <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     ├─GlcNAc (b1-4)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H6N5a6 <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ ├─GlcNAc (b1-6)
  #     │ │ └─Gal (b1-4)
  #     │ └─GlcNAc (b1-2)
  #     │   └─Gal (b1-4)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  make_glycan("Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H6N5S1a6 <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ ├─GlcNAc (b1-6)
  #     │ │ └─Gal (b1-4)
  #     │ └─GlcNAc (b1-2)
  #     │   └─Gal (b1-4)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  #           └─Neu5Ac (a2-3)
  make_glycan("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H6N5a3 <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     │   └─Gal (b1-4)
  #     └─Man (a1-3)
  #       ├─GlcNAc (b1-4)
  #       │ └─Gal (b1-4)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  make_glycan("Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H7N6 <- function(mono_type, linkage) {
  # GlcNAc (b1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ ├─GlcNAc (b1-6)
  #     │ │ └─Gal (b1-4)
  #     │ └─GlcNAc (b1-2)
  #     │   └─Gal (b1-4)
  #     └─Man (a1-3)
  #       ├─GlcNAc (b1-4)
  #       │ └─Gal (b1-4)
  #       └─GlcNAc (b1-2)
  #         └─Gal (b1-4)
  make_glycan("Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", mono_type, linkage)
}

complex_H3N4F1_coreF <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # ├─Fuc (a1-6)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-", mono_type, linkage)
}

complex_H3N4F1_armF <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Fuc (a1-3)
  make_glycan("Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H3N4F2_2coreF <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # ├─Fuc (a1-8)
  # ├─Fuc (a1-6)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  make_glycan("GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)][Fuc(a1-8)]GlcNAc(?1-", mono_type, linkage)
}

complex_H3N4F2_2armF <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     │   └─Fuc (a1-3)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Fuc (a1-3)
  make_glycan("Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-", mono_type, linkage)
}

complex_H3N4F2_1coreF_1armF <- function(mono_type, linkage) {
  # GlcNAc (?1-)
  # ├─Fuc (a1-6)
  # └─GlcNAc (b1-4)
  #   └─Man (b1-4)
  #     ├─Man (a1-6)
  #     │ └─GlcNAc (b1-2)
  #     └─Man (a1-3)
  #       └─GlcNAc (b1-2)
  #         └─Fuc (a1-3)
  make_glycan("Fuc(a1-3)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-", mono_type, linkage)
}

o_glycan_core_1 <- function(mono_type, linkage) {
  make_glycan("Gal(b1-3)GalNAc(a1-", mono_type, linkage)
}


param_grid <- function() {
  expand.grid(
    mono_type = c("generic", "concrete"),
    linkage = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}


# ========== Is N-glycan ==========
patrick::with_parameters_test_that("H4N4 is N-glycan", {
  glycan <- complex_H4N4(mono_type, linkage)
  expect_true(is_n_glycan(glycan))
}, param_grid())


patrick::with_parameters_test_that("O-glycan core 1 is not N-glycan", {
  glycan <- o_glycan_core_1(mono_type, linkage)
  expect_false(is_n_glycan(glycan))
}, param_grid())


test_that("is_n_glycan works vectorizedly", {
  glycans <- c(
    complex_H4N4("generic", linkage = FALSE),
    o_glycan_core_1("generic", linkage = FALSE)
  )
  expect_equal(is_n_glycan(glycans), c(TRUE, FALSE))
})


# ========== N-glycan types ==========
test_that("paucimannose H3N2", {
  glycan <- paucimannose_H3N2("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "paucimannose")
})


test_that("paucimannose H3N2 strict", {
  glycan <- paucimannose_H3N2("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "paucimannose")
})


patrick::with_parameters_test_that("paucimannose H4N2a3", {
  glycan <- paucimannose_H4N2a3(mono_type, linkage)
  expect_identical(n_glycan_type(glycan), "paucimannose")
}, param_grid())


test_that("paucimannose H4N2a3 strict", {
  glycan <- paucimannose_H4N2a3("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "paucimannose")
})


test_that("paucimannose H4N2a6", {
  glycan <- paucimannose_H4N2a6("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "paucimannose")
})


test_that("paucimannose H4N2a6 strict", {
  glycan <- paucimannose_H4N2a6("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "paucimannose")
})


patrick::with_parameters_test_that("highmannose H5N2", {
  glycan <- highmannose_H5N2(mono_type, linkage)
  expect_identical(n_glycan_type(glycan), "highmannose")
}, param_grid())


test_that("highmannose H5N2 strict", {
  glycan <- highmannose_H5N2("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "highmannose")
})


test_that("highmannose H6N2", {
  glycan <- highmannose_H6N2("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "highmannose")
})


test_that("highmannose H6N2 strict", {
  glycan <- highmannose_H6N2("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "highmannose")
})


test_that("hybrid H5N3", {
  glycan <- hybrid_H5N3("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "hybrid")
})


test_that("hybrid H5N3 strict", {
  glycan <- hybrid_H5N3("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "hybrid")
})


patrick::with_parameters_test_that("hybrid H4N3a3", {
  glycan <- hybrid_H4N3a3(mono_type, linkage)
  expect_identical(n_glycan_type(glycan), "hybrid")
}, param_grid())


test_that("hybrid H4N3a3 strict", {
  glycan <- hybrid_H4N3a3("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "hybrid")
})


test_that("hybrid H4N3a3F1", {
  glycan <- hybrid_H4N3a3F1("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "hybrid")
})


test_that("hybrid H4N3a3F1 strict", {
  glycan <- hybrid_H4N3a3F1("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "hybrid")
})


patrick::with_parameters_test_that("complex H3N3", {
  glycan <- complex_H3N3(mono_type, linkage)
  expect_identical(n_glycan_type(glycan), "complex")
}, param_grid())


test_that("complex H3N3 strict", {
  glycan <- complex_H3N3("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "complex")
})


test_that("complex H3N4", {
  glycan <- complex_H3N4("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "complex")
})


test_that("complex H3N4 strict", {
  glycan <- complex_H3N4("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "complex")
})


test_that("complex H4N4", {
  glycan <- complex_H4N4("generic", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "complex")
})


test_that("complex H4N4 strict", {
  glycan <- complex_H4N4("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "complex")
})


test_that("n_glycan_type works vectorizedly", {
  glycans <- c(
    paucimannose_H3N2("generic", linkage = FALSE),
    complex_H4N4("generic", linkage = FALSE)
  )
  expect_equal(n_glycan_type(glycans), c("paucimannose", "complex"))
})


# ========== Bisecting N-glycan ==========
patrick::with_parameters_test_that("complex H3N5 bisect", {
  glycan <- complex_H3N5_bisect(mono_type, linkage)
  expect_true(has_bisecting(glycan))
}, param_grid())


test_that("complex H3N5 bisect strict", {
  glycan <- complex_H3N5_bisect("concrete", linkage = TRUE)
  expect_true(has_bisecting(glycan, strict = TRUE))
})


test_that("complex H4H4 not bisect", {
  glycan <- complex_H4N4("generic", linkage = FALSE)
  expect_false(has_bisecting(glycan))
})


test_that("has_bisecting works vectorizedly", {
  glycans <- c(
    complex_H3N5_bisect("generic", linkage = FALSE),
    complex_H4N4("generic", linkage = FALSE)
  )
  expect_equal(has_bisecting(glycans), c(TRUE, FALSE))
})


# ========== Number of antennae ==========
test_that("one antenna: H3N3", {
  glycan <- complex_H3N3("generic", linkage = FALSE)
  expect_identical(n_antennae(glycan), 1L)
})


patrick::with_parameters_test_that("two antennae: H4N4", {
  glycan <- complex_H4N4(mono_type, linkage)
  expect_identical(n_antennae(glycan), 2L)
}, param_grid())


test_that("two antennae: H3N5, bisecting", {
  glycan <- complex_H3N5_bisect("generic", linkage = FALSE)
  expect_identical(n_antennae(glycan), 2L)
})


test_that("three antennae: H6N5, on a3 mannose", {
  glycan <- complex_H6N5a3("generic", linkage = FALSE)
  expect_identical(n_antennae(glycan), 3L)
})


test_that("three antennae: H6N5, on a3 mannose, strict", {
  glycan <- complex_H6N5a3("concrete", linkage = TRUE)
  expect_identical(n_antennae(glycan, strict = TRUE), 3L)
})


test_that("three antennae: H6N5, on a6 mannose, strict", {
  glycan <- complex_H6N5a6("concrete", linkage = TRUE)
  expect_identical(n_antennae(glycan, strict = TRUE), 3L)
})


test_that("hybrid H5N3 has NA antennae", {
  glycan <- hybrid_H5N3("generic", linkage = FALSE)
  expect_identical(n_antennae(glycan), NA_integer_)
})


test_that("highmannose H5N2 has NA antennae", {
  glycan <- highmannose_H5N2("generic", linkage = FALSE)
  expect_identical(n_antennae(glycan), NA_integer_)
})


test_that("paucimannose H3N2 has NA antennae", {
  glycan <- paucimannose_H3N2("generic", linkage = FALSE)
  expect_identical(n_antennae(glycan), NA_integer_)
})


test_that("n_antennae works vectorizedly", {
  glycans <- c(
    complex_H3N3("generic", linkage = FALSE),
    complex_H4N4("generic", linkage = FALSE)
  )
  expect_equal(n_antennae(glycans), c(1L, 2L))
})


# ========== Number of core fucoses ==========
test_that("no core fucose: H3N3", {
  glycan <- complex_H3N3("generic", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 0L)
})


patrick::with_parameters_test_that("one core fucose: H3N4F1", {
  glycan <- complex_H3N4F1_coreF(mono_type, linkage)
  expect_identical(n_core_fuc(glycan), 1L)
}, param_grid())


test_that("two core fucoses: H3N4F2", {
  glycan <- complex_H3N4F2_2coreF("generic", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 2L)
})


test_that("one core fucose: H3N4F2", {
  glycan <- complex_H3N4F2_1coreF_1armF("generic", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 1L)
})


test_that("one core fucose: H3N4F1, strict", {
  glycan <- complex_H3N4F1_coreF("concrete", linkage = TRUE)
  expect_identical(n_core_fuc(glycan, strict = TRUE), 1L)
})


test_that("no core fucose: H3N4F1", {
  glycan <- complex_H3N4F1_armF("generic", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 0L)
})


test_that("n_core_fuc works vectorizedly", {
  glycans <- c(
    complex_H3N3("generic", linkage = FALSE),
    complex_H3N4F1_coreF("generic", linkage = FALSE)
  )
  expect_equal(n_core_fuc(glycans), c(0L, 1L))
})


# ========== Number of arm fucoses ==========
test_that("no arm fucose: H3N3", {
  glycan <- complex_H3N3("generic", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 0L)
})


patrick::with_parameters_test_that("one arm focuse: H3N4F1", {
  glycan <- complex_H3N4F1_armF(mono_type, linkage)
  expect_identical(n_arm_fuc(glycan), 1L)
}, param_grid())


test_that("two arm fucoses: H3N4F2", {
  glycan <- complex_H3N4F2_2armF("generic", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 2L)
})


test_that("one arm fucose: H3N4F2", {
  glycan <- complex_H3N4F2_1coreF_1armF("generic", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 1L)
})


test_that("one arm fucose: H3N4F1, strict", {
  glycan <- complex_H3N4F1_armF("concrete", linkage = TRUE)
  expect_identical(n_arm_fuc(glycan, strict = TRUE), 1L)
})


test_that("no arm fucose: H3N4F1", {
  glycan <- complex_H3N4F1_coreF("generic", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 0L)
})


test_that("n_arm_fuc works vectorizedly", {
  glycans <- c(
    complex_H3N3("generic", linkage = FALSE),
    complex_H3N4F1_armF("generic", linkage = FALSE)
  )
  expect_equal(n_arm_fuc(glycans), c(0L, 1L))
})


# ========== Number of galactoses ==========
patrick::with_parameters_test_that("one galactose: H4N4S1", {
  glycan <- complex_H4N4S1(mono_type, linkage)
  expect_identical(n_gal(glycan), 1L)
})


test_that("one galactose: H4N4S1, strict", {
  glycan <- complex_H4N4S1("concrete", linkage = TRUE)
  expect_identical(n_gal(glycan, strict = TRUE), 1L)
})


test_that("no galactose: hybrid ", {
  glycan <- highmannose_H5N2("generic", linkage = TRUE)
  expect_identical(n_gal(glycan), 0L)
})


test_that("n_gal works vectorizedly", {
  glycans <- c(
    complex_H4N4S1("generic", linkage = FALSE),
    highmannose_H5N2("generic", linkage = FALSE)
  )
  expect_equal(n_gal(glycans), c(1L, 0L))
})


# ========== Number of terminal galactoses ==========
patrick::with_parameters_test_that("one terminal galactose: H4N4", {
  glycan <- complex_H4N4(mono_type, linkage)
  expect_identical(n_terminal_gal(glycan), 1L)
}, param_grid())


test_that("two terminal galactoses: H6N5S1", {
  glycan <- complex_H6N5S1a6("generic", linkage = FALSE)
  expect_identical(n_terminal_gal(glycan), 2L)
})


test_that("two terminal galactoses: H6N5S1, strict", {
  glycan <- complex_H6N5S1a6("concrete", linkage = TRUE)
  expect_identical(n_terminal_gal(glycan), 2L)
})


test_that("no terminal galactose: H4N4S1", {
  glycan <- complex_H4N4S1("generic", linkage = FALSE)
  expect_identical(n_terminal_gal(glycan), 0L)
})


test_that("no terminal galactose on high-mannose glycan", {
  glycan <- highmannose_H5N2("generic", linkage = FALSE)
  expect_identical(n_terminal_gal(glycan), 0L)
})


test_that("no terminal galactose: hybrid H5N3", {
  glycan <- hybrid_H5N3("generic", linkage = FALSE)
  expect_identical(n_terminal_gal(glycan), 0L)
})


test_that("one terminal galactose: hybrid H6N3", {
  glycan <- hybrid_H6N3("generic", linkage = FALSE)
  expect_identical(n_terminal_gal(glycan), 1L)
})


test_that("n_terminal_gal works vectorizedly", {
  glycans <- c(
    complex_H4N4("generic", linkage = FALSE),
    complex_H6N5S1a6("generic", linkage = FALSE)
  )
  expect_equal(n_terminal_gal(glycans), c(1L, 2L))
})


# ========== Describe N-glycan ==========
test_that("describing glycans works", {
  glycans <- c(
    paucimannose_H3N2("generic", linkage = FALSE),
    highmannose_H5N2("generic", linkage = FALSE),
    complex_H4N4("generic", linkage = FALSE),
    complex_H3N4F2_1coreF_1armF("generic", linkage = FALSE)
  )
  result <- describe_n_glycans(glycans)

  expect_equal(result$glycan_type, c("paucimannose", "highmannose", "complex", "complex"))
  expect_equal(result$bisecting, c(FALSE, FALSE, FALSE, FALSE))
  expect_equal(result$n_antennae, c(NA, NA, 2L, 2L))
  expect_equal(result$n_core_fuc, c(0L, 0L, 0L, 1L))
  expect_equal(result$n_arm_fuc, c(0L, 0L, 0L, 1L))
  expect_equal(result$n_gal, c(0L, 0L, 1L, 0L))
  expect_equal(result$n_terminal_gal, c(0L, 0L, 1L, 0L))
})


test_that("strictly describing glycans works", {
  glycans <- c(
    highmannose_H5N2("concrete", linkage = TRUE),
    complex_H4N4("concrete", linkage = TRUE),
    complex_H3N4F2_1coreF_1armF("concrete", linkage = TRUE)
  )
  expect_snapshot(describe_n_glycans(glycans, strict = TRUE))
})


test_that("describe not N-glycans", {
  glycans <- c(
    o_glycan_core_1("generic", linkage = FALSE),
    highmannose_H5N2("generic", linkage = FALSE)
  )
  expect_error(describe_n_glycans(glycans), "Glycans at indices .* are not N-glycans")
})


test_that("describe N-glycans with no name", {
  glycans <- c(
    highmannose_H5N2("generic", linkage = FALSE),
    complex_H4N4("generic", linkage = FALSE),
    complex_H3N4F2_1coreF_1armF("generic", linkage = FALSE)
  )
  result <- describe_n_glycans(glycans)
  expect_snapshot(result)
})


# ========== Multi-format Support ==========
test_that("N-glycan functions support multiple structure formats", {
  # Test IUPAC-condensed format
  glycan_iupac <- "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  expect_true(is_n_glycan(glycan_iupac))
  expect_equal(n_glycan_type(glycan_iupac), "paucimannose")

  # Test IUPAC-short format
  glycan_short <- "Mana3(Mana6)Manb4GlcNAcb4GlcNAcb-"
  expect_true(is_n_glycan(glycan_short))
  expect_equal(n_glycan_type(glycan_short), "paucimannose")
})


# ========== Vectorized ==========
test_that("vectorized N-glycan properties work", {
  glycans <- c(
    paucimannose_H3N2("generic", linkage = FALSE),
    highmannose_H5N2("generic", linkage = FALSE),
    complex_H4N4("generic", linkage = FALSE)
  )
  expect_equal(n_glycan_type(glycans), c("paucimannose", "highmannose", "complex"))
})

test_that("O-glycans raise warnings", {
  glycans <- c(
    o_glycan_core_1("generic", linkage = FALSE),
    highmannose_H5N2("generic", linkage = FALSE)
  )
  expect_warning(n_glycan_type(glycans), "must be N-glycans")
})