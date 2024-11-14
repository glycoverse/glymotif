# ========== Glycan Graphs ==========
make_glycan <- function(iupac, mono_type, linkage) {
  glycan <- glyparse::parse_iupac_condensed(iupac)
  glycan <- glyrepr::convert_glycan_mono_type(glycan, to = mono_type, strict = FALSE)
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


# ========== N-glycan types ==========
param_grid <- function() {
  expand.grid(
    mono_type = c("simple", "generic", "concrete"),
    linkage = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
}


test_that("paucimannose H3N2", {
  glycan <- paucimannose_H3N2("simple", linkage = FALSE)
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
  glycan <- paucimannose_H4N2a6("simple", linkage = FALSE)
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
  glycan <- highmannose_H6N2("simple", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "highmannose")
})


test_that("highmannose H6N2 strict", {
  glycan <- highmannose_H6N2("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "highmannose")
})


test_that("hybrid H5N3", {
  glycan <- hybrid_H5N3("simple", linkage = FALSE)
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
  glycan <- hybrid_H4N3a3F1("simple", linkage = FALSE)
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
  glycan <- complex_H3N4("simple", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "complex")
})


test_that("complex H3N4", {
  glycan <- complex_H3N4("simple", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "complex")
})


test_that("complex H4N4", {
  glycan <- complex_H4N4("simple", linkage = FALSE)
  expect_identical(n_glycan_type(glycan), "complex")
})


test_that("complex H4N4 strict", {
  glycan <- complex_H4N4("concrete", linkage = TRUE)
  expect_identical(n_glycan_type(glycan, strict = TRUE), "complex")
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
  glycan <- complex_H4N4("simple", linkage = FALSE)
  expect_false(has_bisecting(glycan))
})


# ========== Number of antennae ==========
test_that("one antenna: H3N3", {
  glycan <- complex_H3N3("simple", linkage = FALSE)
  expect_identical(n_antennae(glycan), 1L)
})


patrick::with_parameters_test_that("two antennae: H4N4", {
  glycan <- complex_H4N4(mono_type, linkage)
  expect_identical(n_antennae(glycan), 2L)
}, param_grid())


test_that("two antennae: H3N5, bisecting", {
  glycan <- complex_H3N5_bisect("simple", linkage = FALSE)
  expect_identical(n_antennae(glycan), 2L)
})


test_that("three antennae: H6N5, on a3 mannose", {
  glycan <- complex_H6N5a3("simple", linkage = FALSE)
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
  glycan <- hybrid_H5N3("simple", linkage = FALSE)
  expect_identical(n_antennae(glycan), NA_integer_)
})


test_that("highmannose H5N2 has NA antennae", {
  glycan <- highmannose_H5N2("simple", linkage = FALSE)
  expect_identical(n_antennae(glycan), NA_integer_)
})


test_that("paucimannose H3N2 has NA antennae", {
  glycan <- paucimannose_H3N2("simple", linkage = FALSE)
  expect_identical(n_antennae(glycan), NA_integer_)
})


# ========== Number of core fucoses ==========
test_that("no core fucose: H3N3", {
  glycan <- complex_H3N3("simple", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 0L)
})


patrick::with_parameters_test_that("one core fucose: H3N4F1", {
  glycan <- complex_H3N4F1_coreF(mono_type, linkage)
  expect_identical(n_core_fuc(glycan), 1L)
}, param_grid())


test_that("two core fucoses: H3N4F2", {
  glycan <- complex_H3N4F2_2coreF("simple", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 2L)
})


test_that("one core fucose: H3N4F2", {
  glycan <- complex_H3N4F2_1coreF_1armF("simple", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 1L)
})


test_that("one core fucose: H3N4F1, strict", {
  glycan <- complex_H3N4F1_coreF("concrete", linkage = TRUE)
  expect_identical(n_core_fuc(glycan, strict = TRUE), 1L)
})


test_that("no core fucose: H3N4F1", {
  glycan <- complex_H3N4F1_armF("simple", linkage = FALSE)
  expect_identical(n_core_fuc(glycan), 0L)
})


# ========== Number of arm fucoses ==========
test_that("no arm fucose: H3N3", {
  glycan <- complex_H3N3("simple", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 0L)
})


patrick::with_parameters_test_that("one arm focuse: H3N4F1", {
  glycan <- complex_H3N4F1_armF(mono_type, linkage)
  expect_identical(n_arm_fuc(glycan), 1L)
}, param_grid())


test_that("two arm fucoses: H3N4F2", {
  glycan <- complex_H3N4F2_2armF("simple", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 2L)
})


test_that("one arm fucose: H3N4F2", {
  glycan <- complex_H3N4F2_1coreF_1armF("simple", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 1L)
})


test_that("one arm fucose: H3N4F1, strict", {
  glycan <- complex_H3N4F1_armF("concrete", linkage = TRUE)
  expect_identical(n_arm_fuc(glycan, strict = TRUE), 1L)
})


test_that("no arm fucose: H3N4F1", {
  glycan <- complex_H3N4F1_coreF("simple", linkage = FALSE)
  expect_identical(n_arm_fuc(glycan), 0L)
})
