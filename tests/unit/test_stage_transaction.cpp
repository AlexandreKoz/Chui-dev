// Stage-transaction substrate: reserve -> prepare -> validate -> commit/reject.
// STATUS: written as source only, NOT RUN by the authoring agent (user-reserved).
//
// Invariants:
//  * a governor refusal is recorded (not thrown) so the rank still votes;
//  * prepare() classifies numerical retry vs model downgrade vs resource vs
//    invalid state;
//  * a rejected transaction never runs its commit and leaves the
//    "authoritative" value untouched; reservations are released exactly once;
//  * a peer rejection (collective vote) rejects a locally valid rank;
//  * exactly one vote per transaction; commit before validate is a logic error;
//  * a committed transaction commits pending reservations before the
//    authoritative mutation (accounting matches residency).
#include <stdexcept>
#include <string>

#include "cosmosim/core/memory_governor.hpp"
#include "cosmosim/core/stage_transaction.hpp"
#include "../support/test_require.hpp"

namespace {

namespace co = cosmosim::core;

class FixedPeerDecision final : public co::StageCollectiveDecision {
 public:
  explicit FixedPeerDecision(bool peers_accept) : m_peers_accept(peers_accept) {}
  [[nodiscard]] bool allAccept(bool local_accept, std::string_view) const override {
    ++m_votes;
    return local_accept && m_peers_accept;
  }
  mutable int m_votes = 0;

 private:
  bool m_peers_accept = true;
};

void testCommitPathAndAccounting() {
  co::MemoryGovernor governor(co::MemoryGovernorPolicy{.hard_limit_bytes = 1000U});
  double authoritative = 1.0;
  double staged = 0.0;
  {
    co::StageTransaction transaction("test.commit", &governor);
    CHUI_REQUIRE(transaction.reserve(co::MemoryClass::kScratchArena, 300U, "test.commit.scratch"));
    CHUI_REQUIRE(transaction.prepare([&] { staged = authoritative + 1.0; }));
    CHUI_REQUIRE(authoritative == 1.0);  // prepare never writes truth
    CHUI_REQUIRE(transaction.validate(true, co::StageRejectionKind::kNone, ""));
    const auto outcome = transaction.commit([&] {
      CHUI_REQUIRE(governor.snapshot().committed_bytes == 300U);
      authoritative = staged;
    });
    CHUI_REQUIRE(outcome.committed());
    CHUI_REQUIRE(authoritative == 2.0);
  }
  CHUI_REQUIRE(governor.snapshot().accounted_bytes == 0U);
}

void testResourceRefusalStillVotesAndLeavesTruth() {
  co::MemoryGovernor governor(co::MemoryGovernorPolicy{.hard_limit_bytes = 100U});
  const FixedPeerDecision decision(true);
  double authoritative = 5.0;
  co::StageTransaction transaction("test.refusal", &governor, &decision);
  CHUI_REQUIRE(!transaction.reserve(co::MemoryClass::kScratchArena, 500U, "test.refusal.scratch"));
  CHUI_REQUIRE(!transaction.prepare([&] { authoritative = -1.0; }));  // not executed
  CHUI_REQUIRE(!transaction.validate(true, co::StageRejectionKind::kNone, ""));
  CHUI_REQUIRE(decision.m_votes == 1);
  CHUI_REQUIRE(transaction.outcome().rejection == co::StageRejectionKind::kResourceRefused);
  CHUI_REQUIRE(authoritative == 5.0);
  CHUI_REQUIRE(governor.snapshot().accounted_bytes == 0U);
  bool commit_refused = false;
  try {
    (void)transaction.commit([&] { authoritative = -2.0; });
  } catch (const std::logic_error&) {
    commit_refused = true;
  }
  CHUI_REQUIRE(commit_refused);
  CHUI_REQUIRE(authoritative == 5.0);
}

void testPrepareClassification() {
  {
    co::StageTransaction transaction("test.retry", nullptr);
    CHUI_REQUIRE(!transaction.prepare([] { throw co::StageNumericalRetryRequest("CFL violated"); }));
    CHUI_REQUIRE(transaction.outcome().rejection == co::StageRejectionKind::kNumericalRetry);
  }
  {
    co::StageTransaction transaction("test.downgrade", nullptr);
    CHUI_REQUIRE(!transaction.prepare([] { throw co::StageModelDowngradeRequest("chemistry interval unfinished"); }));
    CHUI_REQUIRE(transaction.outcome().rejection == co::StageRejectionKind::kModelDowngradeRequired);
  }
  {
    co::StageTransaction transaction("test.invalid", nullptr);
    CHUI_REQUIRE(!transaction.prepare([] { throw std::runtime_error("negative density"); }));
    CHUI_REQUIRE(transaction.outcome().rejection == co::StageRejectionKind::kInvalidState);
    CHUI_REQUIRE(transaction.outcome().reason.find("negative density") != std::string::npos);
  }
}

void testPeerRejectionAndSingleVote() {
  const FixedPeerDecision peers_reject(false);
  double authoritative = 3.0;
  co::StageTransaction transaction("test.peer", nullptr, &peers_reject);
  CHUI_REQUIRE(transaction.prepare([] {}));
  CHUI_REQUIRE(!transaction.validate(true, co::StageRejectionKind::kNone, ""));
  CHUI_REQUIRE(transaction.outcome().rejection == co::StageRejectionKind::kPeerRejected);
  CHUI_REQUIRE(authoritative == 3.0);
  bool second_vote_refused = false;
  try {
    (void)transaction.validate(true, co::StageRejectionKind::kNone, "");
  } catch (const std::logic_error&) {
    second_vote_refused = true;
  }
  CHUI_REQUIRE(second_vote_refused);
  CHUI_REQUIRE(peers_reject.m_votes == 1);
}

void testLocalValidationFailureCarriesKind() {
  co::StageTransaction transaction("test.local", nullptr);
  CHUI_REQUIRE(transaction.prepare([] {}));
  CHUI_REQUIRE(!transaction.validate(false, co::StageRejectionKind::kNumericalRetry, "non-finite staged state"));
  CHUI_REQUIRE(transaction.outcome().rejection == co::StageRejectionKind::kNumericalRetry);
  CHUI_REQUIRE(transaction.outcome().phase == co::StageTransactionPhase::kRejected);
}

}  // namespace

int main() {
  testCommitPathAndAccounting();
  testResourceRefusalStillVotesAndLeavesTruth();
  testPrepareClassification();
  testPeerRejectionAndSingleVote();
  testLocalValidationFailureCarriesKind();
  return 0;
}
